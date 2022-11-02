#> Compute energy fractions in longwave bands as a function of temperature
"""
fband = radset()
"""
function radset()

    epslw = 0.05
    fband = fill(NaN, 400, 4)
    eps1 = 1.0 - epslw

    for jtemp = 200: 320
        fband[jtemp, 2] = (0.148 - 3.0e-6 * (jtemp - 247)^2) * eps1
        fband[jtemp, 3] = (0.356 - 5.2e-6 * (jtemp - 282)^2) * eps1
        fband[jtemp, 4] = (0.314 + 1.0e-5 * (jtemp - 315)^2) * eps1
        fband[jtemp, 1] = eps1 - (fband[jtemp, 2] + fband[jtemp, 3] + fband[jtemp, 4])
    end

    for jb = 1: 4
        for jtemp = 100: 199
            fband[jtemp, jb] = fband[200, jb]
        end
        for jtemp = 321: 400
            fband[jtemp, jb] = fband[320, jb]
        end
    end
    return fband
end

"""
ngp = 1
epslw = 0.05
emisfc = 0.98
sbc = 5.67e-8
ta = fill(300., ngp, nlev)
wvi = fill(0.5, nlev, 2)
tau2 = fill(0.5, ngp, nlev, 4)
fband = radset()
fsfcd, dfabs, flux, st4a = radlw_down(ta, fband, epslw, emisfc, sbc, wvi, tau2, ngp)
"""
#> Compute the downward flux of long - wave radiation
function radlw_down(ta::Matrix{Float64}, fband::Matrix{Float64},
    epslw::Float64, emisfc::Float64, sbc::Float64,
    wvi::Matrix{Float64}, tau2::Array{Float64, 3},
    ngp::Int)

    # Parameters
    nlev = 8
    nband = 4
    
    # Derived Parameters
    nl1 = nlev - 1

    # Local
    dfabs = fill(NaN, ngp, nlev)
    fsfcd = fill(NaN, ngp)
    st4a = fill(NaN, ngp, nlev, 2)
    flux = fill(NaN, ngp, 4)

    # 1. Blackbody emission from atmospheric levels.
    #    The linearized gradient of the blakbody emission is computed
    #    from temperatures at layer boundaries, which are interpolated 
    #    assuming a linear dependence of T on log_sigma.
    #    Above the first (top) level, the atmosphere is assumed isothermal.
    # 
    # Temperature at level boundaries
    for k=1:nl1
        for j=1:ngp
            st4a[j, k, 1] = ta[j, k] + wvi[k, 2] * (ta[j, k + 1] - ta[j, k])
        end
    end
    
    # Mean temperature in stratospheric layers
    for j=1:ngp
        st4a[j, 1, 2] = 0.75 * ta[j, 1] + 0.25 * st4a[j, 1, 1]
        st4a[j, 2, 2] = 0.50 * ta[j, 2] + 0.25 * (st4a[j, 1, 1] + st4a[j, 2, 1])
    end

    # Temperature gradient in tropospheric layers
    anis  = 1.0
    anish=0.5*anis
    for k = 3: nl1
        for j=1:ngp
            st4a[j, k, 2] = anish * max(st4a[j, k, 1] - st4a[j, k - 1, 1], 0.)
        end
    end
    for j=1:ngp
        st4a[j, nlev, 2] = anis * max(ta[j, nlev] - st4a[j, nl1, 1], 0.)
    end

    # Blackbody emission in the stratosphere
    for k = 1: 2
        for j=1:ngp
            st4a[j, k, 1] = sbc * st4a[j, k, 2]^4.0
            st4a[j, k, 2] = 0.
        end
    end

    # Blackbody emission in the troposphere
    for k = 3: nlev
        for j=1:ngp
            st3a = sbc * ta[j, k]^3.
            st4a[j, k, 1] = st3a * ta[j, k]
            st4a[j, k, 2] = 4. * st3a * st4a[j, k, 2]
        end
    end

    # 2.0 Initialization of fluxes
    fsfcd .= 0.
    dfabs .= 0.

    # 3.0 Emission ad absorption of longwave downward flux.
    #    For downward emission, a correction term depending on the
    #    local temperature gradient and on the layer transmissivity is
    #    added to the average (full - level) emission of each layer.

    # 3.1  Stratosphere
    k = 1
    for jb = 1: 2
        for j = 1: ngp
            emis = 1. - tau2[j, k, jb]
            brad = fband[convert(Int, ta[j, k]), jb] * (st4a[j, k, 1] + emis * st4a[j, k, 2])
            flux[j, jb] = emis * brad
            dfabs[j, k] = dfabs[j, k] - flux[j, jb]
        end
    end
    flux[:, 3:nband] .= 0.

    # 3.2  Troposphere
    for jb = 1: nband
        for k = 2: nlev
            for j = 1: ngp
                emis = 1. - tau2[j, k, jb]
                brad = fband[convert(Int, ta[j, k]), jb] * (st4a[j, k, 1] + emis * st4a[j, k, 2])
                dfabs[j, k] = dfabs[j, k] + flux[j, jb]
                flux[j, jb] = tau2[j, k, jb] * flux[j, jb] + emis * brad
                dfabs[j, k] = dfabs[j, k] - flux[j, jb]
            end
        end
    end

    # 3.3 Surface downward flux
    for jb = 1: nband
        for j=1:ngp
            fsfcd[j] = fsfcd[j] + emisfc * flux[j, jb]
        end
    end
    eps1 = epslw * emisfc
    for j=1:ngp
        # 3.4 Correction for "black" band (incl. surface reflection)
        corlw = eps1 * st4a[j, nlev, 1]
        dfabs[j, nlev] = dfabs[j, nlev] - corlw
        fsfcd[j] = fsfcd[j] + corlw
    end
    return fsfcd, dfabs, flux, st4a
end

"""
ngp = 1
ts = fill(320., ngp)
fsfcu = compute_bbe(ts, ngp)
"""
function compute_bbe(ts::Vector{Float64}, ngp::Int)
    # Black-body (or grey-body) emission 
    
    emisfc = 0.98
    sbc::Float64 = 5.67e-8

    fsfcu = fill(NaN, ngp)
    esbc=emisfc*sbc
    for j=1:ngp
        tsq=ts[j]*ts[j]
        fsfcu[j]=esbc*tsq*tsq
    end

    return fsfcu
end

"""
ngp = 1
nlev = 8
emisfc = 0.98
epslw = 0.05
ta = fill(300., ngp, nlev)
ts = fill(320., ngp)
fband = radset()
dsig = fill(0.5, nlev)
tau2 = fill(0.5, ngp, nlev, 4)
stratc = fill(0.5, ngp, 2)
fsfcd, dfabs, flux, st4a = radlw_down(ta, fband, epslw, emisfc, sbc, wvi, tau2, ngp)
fsfcu = compute_bbe(ts, ngp)
fsfc, ftop = radlw_up(ta, ts, fsfcd, fsfcu, dfabs, fband, flux, st4a, emisfc, epslw, dsig, tau2, stratc, ngp)
"""
function radlw_up(ta::Matrix{Float64}, ts::Vector{Float64},
    fsfcd::Vector{Float64}, fsfcu::Vector{Float64}, dfabs::Matrix{Float64}, fband::Matrix{Float64},
    flux::Matrix{Float64}, st4a::Array{Float64, 3}, emisfc::Float64, epslw::Float64, dsig::Vector{Float64},
    tau2::Array{Float64, 3}, stratc::Matrix{Float64}, ngp::Int)

    # Parameters
    nlev = 8
    nband = 4
   
    # Local
    fsfc = fill(NaN, ngp)
    ftop = fill(NaN, ngp)

    for j=1:ngp
        fsfc[j] = fsfcu[j] - fsfcd[j]
    end
    
    refsfc = 1.0 - emisfc
    for jb = 1: nband
        for j = 1: ngp
            flux[j, jb] = fband[convert(Int, ts[j]), jb] * fsfcu[j] + refsfc * flux[j, jb]
        end
    end

    # 4.2  Troposphere

    # Correction for "black" band
    for j=1:ngp
        dfabs[j, nlev] = dfabs[j, nlev] + epslw * fsfcu[j]
    end
    for jb = 1: nband
        for k = nlev:-1:2
            for j = 1: ngp
                emis = 1. - tau2[j, k, jb]
                brad = fband[convert(Int, ta[j, k]), jb] * (st4a[j, k, 1] - emis * st4a[j, k, 2])
                dfabs[j, k] = dfabs[j, k] + flux[j, jb]
                flux[j, jb] = tau2[j, k, jb] * flux[j, jb] + emis * brad
                dfabs[j, k] = dfabs[j, k] - flux[j, jb]
            end
        end
    end

    # 4.3  Stratosphere
    k = 1
    for jb = 1: 2
        for j = 1: ngp
            emis = 1.0 - tau2[j, k, jb]
            brad = fband[convert(Int, ta[j, k]), jb] * (st4a[j, k, 1] - emis * st4a[j, k, 2])
            dfabs[j, k] = dfabs[j, k] + flux[j, jb]
            flux[j, jb] = tau2[j, k, jb] * flux[j, jb] + emis * brad
            dfabs[j, k] = dfabs[j, k] - flux[j, jb]
        end
    end

    # Correction for "black" band and polar night cooling
    for j=1:ngp
        corlw1 = dsig[1] * stratc[j, 2] * st4a[j, 1, 1] + stratc[j, 1]
        corlw2 = dsig[2] * stratc[j, 2] * st4a[j, 2, 1]
        dfabs[j, 1] = dfabs[j, 1] - corlw1
        dfabs[j, 2] = dfabs[j, 2] - corlw2
        ftop[j] = corlw1 + corlw2
    end

    # 4.4  Outgoing longwave radiation
    for jb = 1: nband
        for j=1:ngp
            ftop[j] = ftop[j] + flux[j, jb]
        end
    end
    return fsfc, ftop
end


