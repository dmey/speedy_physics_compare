
"""
tyear = 0.5
csol = 6.
lats = [0., 45., 90.]
solar(tyear, csol, lats)
"""
function solar(tyear::Float64, csol::Float64, lats::Vector{Float64})
    # Average daily flux of solar radiation, from Hartmann (1994)

    topsr = fill(NaN, length(lats))
    nlat = length(lats)

    # Compute cosine and sine of latitude
    clat = cos.(lats * pi / 180.0)
    slat = sin.(lats * pi / 180.0)

    # 1.0 Compute declination angle and Earth - Sun distance factor
    pigr = 2.0 * asin(1.0)
    α = 2.0 * pigr * tyear

    ca1 = cos(α)
    sa1 = sin(α)
    ca2 = ca1 * ca1 - sa1 * sa1
    sa2 = 2.0 * sa1 * ca1
    ca3 = ca1 * ca2 - sa1 * sa2
    sa3 = sa1 * ca2 + sa2 * ca1

    decl = 0.006918 - 0.399912 * ca1 + 0.070257 * sa1 - 0.006758 * ca2 + 0.000907 * sa2 - 0.002697 * ca3 + 0.001480 * sa3

    fdis = 1.000110 + 0.034221 * ca1 + 0.001280 * sa1 + 0.000719 * ca2 + 0.000077 * sa2

    cdecl = cos(decl)
    sdecl = sin(decl)
    tdecl = sdecl / cdecl

    # 2.0 Compute daily - average insolation at the atm. top
    csolp = csol / pigr

    for j in 1:nlat
        ch0 = min(1.0, max(-1.0, -tdecl * slat[j] / clat[j]))
        h0 = acos(ch0)
        sh0 = sin(h0)

        topsr[j] = csolp * fdis * (h0 * slat[j] * sdecl + sh0 * clat[j] * cdecl)
    end
    return topsr

end


"""
tyear = 0.5
nlon = 2
lats = [0., 89.]
sol_oz(tyear, nlon, lats)
"""
function sol_oz(tyear::Float64, nlon::Int, lats::Vector{Float64})

    solc = 342.0 ## Solar constant (area averaged) in W / m^2
    epssw = 0.020 ## Fraction of incoming solar radiation absorbed by ozone

    fsol = fill(NaN, nlon * length(lats))
    ozupp = fill(NaN, nlon * length(lats))
    ozone = fill(NaN, nlon * length(lats))
    zenit = fill(NaN, nlon * length(lats))
    stratz = fill(NaN, nlon * length(lats))

    nlat = length(lats)
    ngp = nlon * nlat

    # Compute cosine and sine of latitude
    clat = cos.(lats * π / 180.0)
    slat = sin.(lats * π / 180.0)

    # α = year phase ( 0 - 2pi, 0 = winter solstice = 22dec.h00 )
    α = 4.0 * asin(1.0) * (tyear + 10.0 / 365.0)
    dα = 0.0

    coz1 = 1.0 * max(0.0, cos(α - dα))
    coz2 = 1.8
    azen = 1.0
    nzen = 2

    rzen = -cos(α) * 23.45 * asin(1.0) / 90.0
    czen = cos(rzen)
    szen = sin(rzen)

    fs0 = 6.0

    # # Solar radiation at the top
    topsr = solar(tyear, 4.0 * solc, lats)

    for j in 1:nlat
        j0 = 1 + nlon * (j - 1)
        flat2 = 1.5 * slat[j]^2 - 0.5

        # Solar radiation at the top
        fsol[j0] = topsr[j]

        # Ozone depth in upper and lower stratosphere
        ozupp[j0] = 0.5 * epssw
        ozone[j0] = 0.4 * epssw * (1.0 + coz1 * slat[j] + coz2 * flat2)

        # Zenith angle correction to (downward) absorptivity
        zenit[j0] = 1.0 + azen * (1.0 - (clat[j] * czen + slat[j] * szen))^nzen

        # Ozone absorption in upper and lower stratosphere
        ozupp[j0] = fsol[j0] * ozupp[j0] * zenit[j0]
        ozone[j0] = fsol[j0] * ozone[j0] * zenit[j0]

        # Polar night cooling in the stratosphere
        stratz[j0] = max(fs0 - fsol[j0], 0.0)

        for i in 1:nlon-1
            fsol[i+j0] = fsol[j0]
            ozone[i+j0] = ozone[j0]
            ozupp[i+j0] = ozupp[j0]
            zenit[i+j0] = zenit[j0]
            stratz[i+j0] = stratz[j0]
        end
    end

    return topsr, fsol, ozupp, ozone, zenit, stratz
end


"""
ngp = 2
nlev = 2

prog_q = fill(1., ngp, nlev)
rh = fill([0.5 0.5]; [0.1 0.1], ngp, nlev)
rh = [0.5  0.2
      0.4  0.1]
cnv_prec = fill(1., ngp)
precls = fill(1., ngp)
cnv_top = fill(1, ngp)
gse = fill(1., ngp)
fmask = fill(1., ngp)

cloud(prog_q, rh, cnv_prec, precls, cnv_top, gse, fmask, ngp, nlev)
"""
function cloud(prog_q::Matrix{Float64}, rh::Matrix{Float64}, cnv_prec::Vector{Float64}, precls::Vector{Float64},
    cnv_top::Vector{Int}, gse::Vector{Float64}, fmask::Vector{Float64}, ngp::Int, nlev::Int)

    rhcl1 = 0.30 # Relative humidity threshold corresponding to cloud cover = 0
    rhcl2 = 1.00 # Relative humidity correponding to cloud cover = 1
    qcl = 0.20 # Specific humidity threshold for cloud cover
    pmaxcl = 10.0 # Maximum value of precipitation (mm/day) contributing to cloud cover
    wpcl = 0.2  # Cloud cover weight for the square-root of precipitation (for p = 1 mm/day)
    gse_s1 = 0.40 # Gradient of dry static energy corresponding to stratiform cloud cover = 1
    gse_s0 = 0.25 # Gradient of dry static energy corresponding to stratiform cloud cover = 0
    clsmax = 0.60 # Maximum stratiform cloud cover
    clsminl = 0.15 # Minimum stratiform cloud cover over land (for RH = 1)

    nl1 = nlev - 1
    nlp = nlev + 1
    rrcl = 1.0 / (rhcl2 - rhcl1)

    cloudc = fill(NaN, ngp)
    icltop = fill(NaN, ngp)
    qcloud = fill(NaN, ngp)
    clstr = fill(NaN, ngp)

    # 1.0 Cloud cover, defined as the sum of:
    # - a term proportional to the square - root of precip. rate
    # - a quadratic function of the max. relative humidity
    #    in tropospheric layers above pbl where q > prog_qcl :
    #     ( = 0 for rhmax < rhcl1, = 1 for rhmax > rhcl2 )
    #    Cloud - top level: defined as the highest (i.e. least sigma)
    #      between the top of convection / condensation and
    #      the level of maximum relative humidity.
    for j in 1:ngp
        if (rh[j, nl1] > rhcl1)
            cloudc[j] = rh[j, nl1] - rhcl1
            icltop[j] = nl1
        else
            cloudc[j] = 0.0
            icltop[j] = nlp
        end
    end

    for k in 3:nlev-2
        for j in 1:ngp
            drh = rh[j, k] - rhcl1
            if (drh > cloudc[j]) & (prog_q[j, k] > qcl)
                cloudc[j] = drh
                icltop[j] = k
            end
        end
    end

    for j in 1:ngp
        pr1 = min(pmaxcl, 86.4 * (cnv_prec[j] + precls[j]))
        cloudc[j] = min(1.0, wpcl * sqrt(pr1) + min(1.0, cloudc[j] * rrcl)^2.0)
        icltop[j] = min(cnv_top[j], icltop[j])
    end

    # 2.0  Equivalent specific humidity of clouds
    for j = 1:ngp
        qcloud[j] = prog_q[j, nl1]
    end

    # 3.0 Stratiform clouds at the top of pbl
    clfact = 1.2
    rgse = 1.0 / (gse_s1 - gse_s0)

    for j in 1:ngp
        # Stratocumulus clouds over sea
        fstab = max(0.0, min(1.0, rgse * (gse[j] - gse_s0)))
        clstr[j] = fstab * max(clsmax - clfact * cloudc[j], 0.0)

        # Stratocumulus clouds over land
        clstrl = max(clstr[j], clsminl) * rh[j, nlev]
        clstr[j] = clstr[j] + fmask[j] * (clstrl - clstr[j])
    end
    return icltop, cloudc, clstr, qcloud

end


"""
ngp = 2
nlev = 3

prog_sp = fill(1., ngp)
prog_q = fill(1., ngp, nlev)
icltop = fill(1, ngp)
cloudc = fill(1., ngp)
clstr = fill(1., ngp)
ozupp = fill(1., ngp)
ozone = fill(1., ngp)
zenit = fill(1., ngp)
stratz = fill(1., ngp)
fsol = fill(1., ngp)
qcloud = fill(1., ngp)
albsfc = fill(1., ngp)
sig = fill(1., nlev)
dsig = fill(1., nlev)

radsw(prog_sp, prog_q, icltop, cloudc, clstr, ozupp, ozone, zenit, stratz, fsol,
    qcloud, albsfc, ngp, nlev, sig, dsig);
"""
function radsw(prog_sp::Vector{Float64}, prog_q::Matrix{Float64}, icltop::Vector{Int}, cloudc::Vector{Float64}, clstr::Vector{Float64},
    ozupp::Vector{Float64}, ozone::Vector{Float64}, zenit::Vector{Float64}, stratz::Vector{Float64}, fsol::Vector{Float64},
    qcloud::Vector{Float64}, albsfc::Vector{Float64}, ngp::Int, nlev::Int, sig::Vector{Float64}, dsig::Vector{Float64})

    albcl = 0.43  # Cloud albedo (for cloud cover = 1)
    albcls = 0.50  # Stratiform cloud albedo (for st. cloud cover = 1)
    abscl1 = 0.015 # Absorptivity of clouds (visible band, maximum value)
    abscl2 = 0.15  # Absorptivity of clouds (visible band, for dq_base = 1 g/kg)

    absdry = 0.033 # Absorptivity of dry air (visible band)
    absaer = 0.033 # Absorptivity of aerosols (visible band)
    abswv1 = 0.022 # Absorptivity of water vapour (visible band, for dq = 1 g/kg)
    abswv2 = 15.0  # Absorptivity of water vapour (near IR band, for dq = 1 g/kg)

    ablwin = 0.3   # Absorptivity of air in "window" band
    ablco2 = 6.0   # Absorptivity of air in CO2 band
    ablwv1 = 0.7   # Absorptivity of water vapour in H2O band 1 (weak), (for dq = 1 g/kg)
    ablwv2 = 50.0   # Absorptivity of water vapour in H2O band 2 (strong), (for dq = 1 g/kg)
    ablcl2 = 0.6   # Absorptivity of "thin" upper clouds in window and H2O bands
    ablcl1 = 12.0   # Absorptivity of "thick" clouds in window band (below cloud top)
    epslw = 0.05    # Fraction of blackbody spectrum absorbed/emitted by PBL only

    ssrd = fill(NaN, ngp)
    ssr = fill(NaN, ngp)
    tsr = fill(NaN, ngp)
    tau2 = fill(NaN, ngp, nlev, 4)
    tend_t_rsw = fill(NaN, ngp, nlev)

    acloud = fill(NaN, ngp)
    psaz = fill(NaN, ngp)
    flux = fill(NaN, ngp, 2)
    stratc = fill(NaN, ngp, 2)

    #
    #      equivalence (frefl(1, 1), tau2(1, 1, 3))

    nl1 = nlev - 1

    fband2 = 0.05
    fband1 = 1.0 - fband2

    # 1.0 Initialization
    tau2 = fill!(tau2, 0)

    for j in 1:ngp
        # Change to ensure only ICLTOP < = NLEV used
        if (icltop[j] <= nlev)
            tau2[j, icltop[j], 3] = albcl * cloudc[j]
        end
        tau2[j, nlev, 3] = albcls * clstr[j]
    end
    # 2.0 Shortwave transmissivity:
    # function of layer mass, ozone (in the statosphere),
    # abs. humidity and cloud cover (in the troposphere)
    for j in 1:ngp
        psaz[j] = prog_sp[j] * zenit[j]
        acloud[j] = cloudc[j] * min(abscl1 * qcloud[j], abscl2)
    end

    for j in 1:ngp
        deltap = psaz[j] * dsig[1]
        tau2[j, 1, 1] = exp(-deltap * absdry)
    end

    for k in 2:nl1
        abs1 = absdry + absaer * sig[k]^2
        for j in 1:ngp
            deltap = psaz[j] * dsig[k]
            if (k >= icltop[j])
                tau2[j, k, 1] = exp(-deltap * (abs1 + abswv1 * prog_q[j, k] + acloud[j]))
            else
                tau2[j, k, 1] = exp(-deltap * (abs1 + abswv1 * prog_q[j, k]))
            end
        end
    end

    abs1 = absdry + absaer * sig[nlev]^2
    for j in 1:ngp
        deltap = psaz[j] * dsig[nlev]
        tau2[j, nlev, 1] = exp(-deltap * (abs1 + abswv1 * prog_q[j, nlev]))
    end

    for k in 2:nlev
        for j in 1:ngp
            deltap = psaz[j] * dsig[k]
            tau2[j, k, 2] = exp(-deltap * abswv2 * prog_q[j, k])
        end
    end

    # 3.0 Shortwave downward flux
    # 3.1 Initialization of fluxes
    for j in 1:ngp
        tsr[j] = fsol[j]
        flux[j, 1] = fsol[j] * fband1
        flux[j, 2] = fsol[j] * fband2
    end

    # 3.2 Ozone and dry - air absorption in the stratosphere
    k = 1
    for j in 1:ngp
        tend_t_rsw[j, k] = flux[j, 1]
        flux[j, 1] = tau2[j, k, 1] * (flux[j, 1] - ozupp[j] * prog_sp[j])
        tend_t_rsw[j, k] = tend_t_rsw[j, k] - flux[j, 1]
    end

    k = 2
    for j = 1:ngp
        tend_t_rsw[j, k] = flux[j, 1]
        flux[j, 1] = tau2[j, k, 1] * (flux[j, 1] - ozone[j] * prog_sp[j])
        tend_t_rsw[j, k] = tend_t_rsw[j, k] - flux[j, 1]
    end

    # 3.3  Absorption and reflection in the troposphere
    for k in 3:nlev
        for j in 1:ngp
            tau2[j, k, 3] = flux[j, 1] * tau2[j, k, 3]
            flux[j, 1] = flux[j, 1] - tau2[j, k, 3]
            tend_t_rsw[j, k] = flux[j, 1]
            flux[j, 1] = tau2[j, k, 1] * flux[j, 1]
            tend_t_rsw[j, k] = tend_t_rsw[j, k] - flux[j, 1]
        end
    end

    for k in 2:nlev
        for j in 1:ngp
            tend_t_rsw[j, k] = tend_t_rsw[j, k] + flux[j, 2]
            flux[j, 2] = tau2[j, k, 2] * flux[j, 2]
            tend_t_rsw[j, k] = tend_t_rsw[j, k] - flux[j, 2]
        end
    end


    # 4.0 Shortwave upward flux
    # 4.1  Absorption and reflection at the surface
    for j in 1:ngp
        ssrd[j] = flux[j, 1] + flux[j, 2]
        flux[j, 1] = flux[j, 1] * albsfc[j]
        ssr[j] = ssrd[j] - flux[j, 1]
    end

    # 4.2  Absorption of upward flux
    for k in nlev:-1:1
        for j in 1:ngp
            tend_t_rsw[j, k] = tend_t_rsw[j, k] + flux[j, 1]
            flux[j, 1] = tau2[j, k, 1] * flux[j, 1]
            tend_t_rsw[j, k] = tend_t_rsw[j, k] - flux[j, 1]
            flux[j, 1] = flux[j, 1] + tau2[j, k, 3]
        end
    end

    # 4.3  Net solar radiation = incoming - outgoing
    for j in 1:ngp
        tsr[j] = tsr[j] - flux[j, 1]
    end

    # 5.0  Initialization of longwave radiation model
    # 5.1  Longwave transmissivity:
    #      function of layer mass, abs. humidity and cloud cover.

    #    Cloud - free levels (stratosphere + PBL)
    k = 1
    for j in 1:ngp
        deltap = prog_sp[j] * dsig[k]
        tau2[j, k, 1] = exp(-deltap * ablwin)
        tau2[j, k, 2] = exp(-deltap * ablco2)
        tau2[j, k, 3] = 1.0
        tau2[j, k, 4] = 1.0
    end

    for k in 2:nlev-2:nlev
        for j in 1:ngp
            deltap = prog_sp[j] * dsig[k]
            tau2[j, k, 1] = exp(-deltap * ablwin)
            tau2[j, k, 2] = exp(-deltap * ablco2)
            tau2[j, k, 3] = exp(-deltap * ablwv1 * prog_q[j, k])
            tau2[j, k, 4] = exp(-deltap * ablwv2 * prog_q[j, k])
        end
    end

    # Cloudy layers (free troposphere)
    for j in 1:ngp
        acloud[j] = cloudc[j] * ablcl2
    end

    for k in 3:nl1
        for j in 1:ngp
            deltap = prog_sp[j] * dsig[k]
            if (k < icltop[j])
                acloud1 = acloud[j]
            else
                acloud1 = ablcl1 * cloudc[j]
            end
            tau2[j, k, 1] = exp(-deltap * (ablwin + acloud1))
            tau2[j, k, 2] = exp(-deltap * ablco2)
            tau2[j, k, 3] = exp(-deltap * max(ablwv1 * prog_q[j, k], acloud[j]))
            tau2[j, k, 4] = exp(-deltap * max(ablwv2 * prog_q[j, k], acloud[j]))
        end
    end

    # 5.2  Stratospheric correction terms
    eps1 = epslw / (dsig[1] + dsig[2])
    for j in 1:ngp
        stratc[j, 1] = stratz[j] * prog_sp[j]
        stratc[j, 2] = eps1 * prog_sp[j]
    end

    return ssrd, ssr, tsr, tau2, tend_t_rsw

end
