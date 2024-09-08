struct Thermo{
   //partition functions
    pftot:   f64,
    pfelec:  f64,
    pftrans: f64,
    pfrot:   f64,
    pfvib:   f64,

   //entropy functions
    stot:    f64,
    selec:   f64,
    strans:  f64,
    srot:    f64,
    svib:    f64,

   //internal energy functions
    utot:    f64,
    uelec:   f64,
    utrans:  f64,
    urot:    f64,
    uvib:    f64,

   //enthalpy functions
    htot:    f64,
    helec:   f64,
    htrans:  f64,
    hrot:    f64,
    hvib:    f64,

   //Helmholtz free energy functions
    ftot:    f64,
    felec:   f64,
    ftrans:  f64,
    frot:    f64,
    fvib:    f64,

   //Gibbs free energy functions
    gtot:    f64,
    gelec:   f64,
    gtrans:  f64,
    grot:    f64,
    gvib:    f64 
}

