"""Reproducible injection/recovery tests; no external observations required.

Run: python examples/recovery_validation.py --output validation/recovery
All presets are fixed before these held-out realizations are generated.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.ndimage import gaussian_filter
from astropy.wcs import WCS
from gtlmapping import GTLMapper
from gtlmapping.exceptions import InsufficientSamplesError

PROFILES=("bt12","conservative","moderate","liberal")
SCENES=("constant","linear","quadratic","no_opaque","unmodeled_structure",
        "background_high","background_low","high_noise","correlated_noise","broad_psf","few_cores")


def scene(seed,kind):
    rng=np.random.default_rng(seed)
    y,x=np.indices((128,128),dtype=float);xn=(x-63.5)/63.5;yn=(y-63.5)/63.5
    tau=.2+.35*np.exp(gaussian_filter(rng.normal(size=x.shape),8))
    tau+=1.3*np.exp(-.5*((x-64-9*np.sin(y/25))/8)**2)
    if kind!="no_opaque":
        locations=np.linspace(17,110,2 if kind=="few_cores" else 4)
        for cy in locations:
            for cx in locations:
                xc,yc=cx+rng.uniform(-4,4),cy+rng.uniform(-4,4)
                width=rng.uniform(3.5,5.5)
                tau+=rng.uniform(7,10)*np.exp(-((x-xc)**2+(y-yc)**2)/(2*width**2))
    fg=np.full(x.shape,20.)
    if kind not in ("constant","no_opaque"):fg+=5*xn+3*yn
    if kind=="quadratic":fg+=3*xn*yn+2*yn**2
    if kind=="unmodeled_structure":fg+=4*np.sin(2*np.pi*x/128)*np.cos(2*np.pi*y/128)
    bg=np.full(x.shape,90.)
    psf=2.5 if kind=="broad_psf" else 1.
    noiseless=gaussian_filter(fg+(bg-fg)*np.exp(-tau),psf,mode="reflect")
    fg=gaussian_filter(fg,psf,mode="reflect")
    trans=noiseless-fg
    true_sigma=-np.log(trans/(bg-fg))/7.5
    noise_rms=1.2 if kind=="high_noise" else .6
    noise=gaussian_filter(rng.normal(size=x.shape),2. if kind=="correlated_noise" else .8,mode="reflect")
    noise*=noise_rms/np.std(noise)
    assumed_bg=bg*(1.1 if kind=="background_high" else .9 if kind=="background_low" else 1.)
    wcs=WCS(naxis=2);wcs.wcs.ctype=["GLON-TAN","GLAT-TAN"]
    wcs.wcs.crval=[28.37,.07];wcs.wcs.crpix=[64.5,64.5];wcs.wcs.cdelt=[-.6/3600,.6/3600]
    return noiseless+noise,assumed_bg,fg,true_sigma,trans,noise_rms,wcs


def run(output,seeds=20):
    output.mkdir(parents=True,exist_ok=True)
    records=[]
    for kind in SCENES:
        for index in range(seeds):
            obs,bg,truth,true_sigma,true_trans,noise,wcs=scene(18100+index,kind)
            mask=np.ones(obs.shape,bool)
            detector=GTLMapper(obs,wcs=wcs)
            samples=detector.detect_foreground(noise_sigma=noise)
            for profile in PROFILES:
                row={"scene":kind,"seed":18100+index,"profile":profile,"samples":len(samples)}
                mapper=GTLMapper(obs,wcs=wcs);mapper.set_background(bg)
                try:
                    fit=mapper.fit_foreground(method=profile,samples=samples,region_mask=mask,
                                              noise_sigma=noise,min_separation_arcsec=8.)
                    # Every profile uses the same adopted image-noise sensitivity.
                    if profile in ("moderate","liberal"):
                        result=getattr(mapper,"compute_"+profile)(intensity_floor=2*noise,bright_pixel_policy="zero")
                    else:
                        result=mapper.compute(intensity_floor=2*noise,bright_pixel_policy="zero")
                    selected=true_trans>3*noise
                    finite=~np.ma.getmaskarray(result.surface_density)
                    if not np.all(finite[selected]):raise ValueError("Missing product on truth-detectable pixels")
                    sigma=result.surface_density.data
                    detected=selected & result.detection_mask
                    row.update(status="ok",foreground_rmse=float(np.sqrt(np.mean((result.foreground-truth)**2))),
                        sum_bias=float(sigma[selected].sum()/true_sigma[selected].sum()-1),
                        sigma_rmse=float(np.sqrt(np.mean((sigma[selected]-true_sigma[selected])**2))),
                        unresolved_fraction=float(result.unresolved_mask.mean()),
                        false_unresolved_fraction=float(np.mean(result.unresolved_mask[selected])),
                        truth_detectable_fraction=float(selected.mean()),
                        detected_fraction_of_truth_detectable=float(detected.sum()/selected.sum()),
                        strict_count=int(result.saturated_mask.sum()),
                        limit_confidence="not_calibrated")
                except (ValueError,RuntimeError,InsufficientSamplesError) as exc:
                    row.update(status="failed",error=f"{type(exc).__name__}: {exc}")
                records.append(row)
        (output/"trials.json").write_text(json.dumps(records,indent=2,allow_nan=False)+"\n")
        print(kind,flush=True)
    summary=[]
    for kind in SCENES:
        for profile in PROFILES:
            rows=[r for r in records if r["scene"]==kind and r["profile"]==profile]
            ok=[r for r in rows if r["status"]=="ok"]
            row={"scene":kind,"profile":profile,"attempted":len(rows),"succeeded":len(ok)}
            for metric in ("foreground_rmse","sum_bias","sigma_rmse","unresolved_fraction","false_unresolved_fraction"):
                if ok:
                    values=[r[metric] for r in ok]
                    for name,q in (("p16",.16),("median",.5),("p84",.84)):
                        row[f"{metric}_{name}"]=float(np.quantile(values,q))
            summary.append(row)
    (output/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path("validation/recovery"))
    parser.add_argument("--seeds",type=int,default=20)
    args=parser.parse_args();run(args.output,args.seeds)
