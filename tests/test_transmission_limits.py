import numpy as np
import pytest
from astropy.io import fits
from gtlmapping import GTLMapper
from gtlmapping.models import InterpolationResult, ForegroundSamples
from gtlmapping.extinction import compute_extinction, unresolved_transmission, transmission_std
from gtlmapping.foreground import _trend_support
from gtlmapping.exceptions import InsufficientSamplesError
from gtlmapping.opacity import FilterOpacity, get_filter_opacity


def test_zero_crossing_and_threshold_are_continuous_in_display_value():
    fg = np.array([28.799999,28.8,28.800001,29.99,30.0,30.01])
    obs=np.full(fg.shape,30.);bg=np.full(fg.shape,70.)
    tau, _, strict, _, _ = compute_extinction(obs,bg,fg,saturation_policy="lower_limit",intensity_floor=1.2)
    limits=unresolved_transmission(obs,fg,1.2)
    assert limits.tolist()==[False,True,True,True,True,True]
    assert strict.tolist()==[False,False,False,False,True,True]
    assert np.max(np.abs(np.diff(tau[:3]))) < 1e-5
    assert np.all(np.isfinite(tau))
    assert np.max(tau) < 4


def test_machine_precision_residual_is_not_a_detection():
    obs=np.nextafter(30.,np.inf)
    tau,_,strict,_,_=compute_extinction(obs,70.,30.,saturation_policy="lower_limit",intensity_floor=1.2)
    assert not strict
    assert unresolved_transmission(obs,30.,1.2)
    assert tau == pytest.approx(-np.log(1.2/40))


def test_mask_policy_excludes_positive_weak_transmission():
    out=compute_extinction([30.,31.,32.],[70.]*3,[30.]*3,detection_threshold=1.2)
    assert out[0].mask.tolist()==[True,True,False]


def test_spatial_threshold_and_invalid_inputs():
    out=compute_extinction([30.,31.,32.,np.nan],[70.,70.,29.,70.],[30.]*4,
                          saturation_policy="lower_limit",intensity_floor=[1.2,np.nan,1.2,1.2])
    assert out[0].mask.tolist()==[False,True,True,True]
    assert out[3][2]
    missing = compute_extinction([31.], [70.], [30.],
                                 saturation_policy="lower_limit", intensity_floor=1.2,
                                 detection_threshold=np.nan)
    assert missing[0].mask.all()


@pytest.mark.parametrize("bad",[0.,-1.,np.inf])
def test_invalid_threshold_rejected(bad):
    with pytest.raises(ValueError): unresolved_transmission([2.],[1.],bad)


def test_mapper_masks_limit_uncertainty_and_writes_separate_masks(tmp_path,galactic_header):
    mapper=GTLMapper(np.array([[30.,30.6,35.]]),header=galactic_header)
    mapper.foreground_result=InterpolationResult(np.full((1,3),30.),"bt12",diagnostics={"noise_sigma":.6})
    mapper.set_background(np.full((1,3),70.))
    result=mapper.compute(observed_std=.6)
    assert result.saturated_mask.tolist()==[[True,False,False]]
    assert result.unresolved_mask.tolist()==[[True,True,False]]
    assert result.detection_mask.tolist()==[[False,False,True]]
    assert result.uncertainty.surface_density_std.mask.tolist()==[[True,True,False]]
    assert not result.surface_density.mask.any()
    with fits.open(result.write(tmp_path/"limits.fits")) as hdus:
        assert np.array_equal(hdus["UNRESOLVED"].data,result.unresolved_mask)
        assert np.array_equal(hdus["SATURATED"].data,result.saturated_mask)
        assert np.allclose(hdus["TRANS_LIM"].data,1.2)


def test_quadratic_requires_full_design_rank():
    x=np.arange(6);y=x*x
    samples=ForegroundSamples(y,x,np.ones(6),np.ones(6,int),6,0,6,1.,1.)
    design=np.column_stack([np.ones(6),x,y,x*x,x*y,y*y])
    with pytest.raises(InsufficientSamplesError,match="design rank 6"):
        _trend_support(design,samples,np.ones((30,6),bool))


def test_residual_covariance_and_mass_conventions():
    assert transmission_std(2.,1.,.5)==pytest.approx(2.)
    with pytest.raises(ValueError): transmission_std(1.,1.,2.)
    assert get_filter_opacity("F480M",mass_basis="total")==pytest.approx(9.76*156/157)
    ref=FilterOpacity("F480M",4.8,"example",9.76,156.,"total")
    assert ref.at_gas_to_dust_ratio(100,mass_basis="total")==pytest.approx(9.76*157/101)
    with pytest.raises(ValueError): get_filter_opacity("F480M",gas_to_dust_ratio=np.nan)
