f.pbk <- function(t, A, parms, ...)
{
  #  Physiologically Based Kinetic Model structure for T-2, HT-2 and their metabolites. 
  #  
  #  written by J. Westerhout, PhD (joost.westerhout@rivm.nl)
  #  National Institute for Public Health and the Environment (RIVM, NL)
  #
  #  Based on PhysB model from Bosgra et al., 2012. DOI: 10.3109/10408444.2012.709225
  
  dots <- list(...)
  
  # Correct for activity level
  if (T) {
  parms[c(2,20:35)] <- lapply(parms[c(2,20:35)], dots$f.select, index=(t/dots$tstep + 1))}
  
  with(as.list(c(A, parms, dots)),
       {
         # Exposed compartments
         
         # Perfusion depends on activity level. Select perfusion at current time point.
        
         if (F) { 
           Qc = Qc(t/tstep + 1)
           Qhrt = Qhrt(t/tstep + 1)
           Qskn = Qskn(t/tstep + 1)
           Qadp = Qadp(t/tstep + 1)
           Qmus = Qmus(t/tstep + 1)
           Qbon = Qbon(t/tstep + 1)
           Qbrn = Qbrn(t/tstep + 1)
           Qthy = Qthy(t/tstep + 1)
           Qgon = Qgon(t/tstep + 1)
           Qkid = Qkid(t/tstep + 1)
           Qsto = Qsto(t/tstep + 1)
           Qint = Qint(t/tstep + 1)
           Qspl = Qspl(t/tstep + 1)
           Qpan = Qpan(t/tstep + 1)
           Qliv = Qliv(t/tstep + 1)
           Qrem = Qrem(t/tstep + 1)
           Qp   = Qp(t/tstep + 1)}
         
         #### T2 ----
         ## Inhalation 
         Cair_T2 <- 0 # X.air_T2[t/tstep + 1]
         Calv_T2 <- Aalv_T2/Valv
         
         Vv <- 0.67*Vb
         Va <- 0.33*Vb
         Cv_T2 <- Av_T2/Vv
         Ca_T2 <- Aa_T2/Va
         
         if (include_lung == TRUE) {
           Clun_T2 <- Alun_T2/Vlun
           dAalv_T2 <- Qp*(Cair_T2 - Calv_T2) - PAalv_T2*(Calv_T2 - Clun_T2/Pluna_T2)     # Alveolar air   
           dAlun_T2 <- Qc*(Cv_T2 - Clun_T2/Plunb_T2) + PAalv_T2*(Calv_T2 - Clun_T2/Pluna_T2) # Lungs
           Aair_T2 <- Qp*(Cair_T2 - Calv_T2) # outside air
         }
         
         else {
           Clun_T2 <- Aalv_T2/Valv
           dAalv_T2 <- Qp*(Cair_T2 - Calv_T2) + Qc*(Cv_T2 - Calv_T2/Pab_T2)  # Alveolar air
           dAlun_T2 <- 0 # excluded
           Aair_T2 <- Qp*(Cair_T2 - Calv_T2) # outside air
         }

         ## Skin contact
         Css_T2 <- Ass_T2/Vss # Vss in L (dm3)
         
         As_T2 <- 0 # X.skin_T2[t/tstep + 1] # amount of substance added to the skin at time t
         dAss_T2  <- As_T2 - Ks_T2*fss*BSA*Css_T2/100                             # Skin surface; divided by 100 to correct for BSA unit (m2 -> dm2)
         
         QskE <- fss*Qskn 
         
         CskE_T2 <- ifelse(fss!=0, AskE_T2/(fss*Vskn), 0)
         dAskE_T2 <- Ks_T2*fss*BSA*Css_T2/100 + QskE*(Ca_T2 - CskE_T2/Pskn_T2)      # Exposed skin
         
         QskU <- (1 - fss)*Qskn
         CskU_T2 <- AskU_T2/((1 - fss)*Vskn)
         dAskU_T2 <- QskU*(Ca_T2 - CskU_T2/Pskn_T2)                       # Unexposed skin
         
         ## Oral
         #Aoral_T2 <- X.oral_T2[t/tstep + 1]
         Aoral_T2 <- X.oral_T2(t) 
         
         dAoral_total_T2 <- Aoral_T2
         dAstoL_T2 <- Aoral_T2 - ka_sto_T2*AstoL_T2 - kse*AstoL_T2 - 
           km_sto_T2_HT2*AstoL_T2 - km_sto_T2_T2_Phase_II*AstoL_T2 - km_sto_T2_T2_Phase_I*AstoL_T2             # Stomach lumen
         dAintL_T2 <- kse*AstoL_T2 - ka_int_T2*AintL_T2 - kfec*AintL_T2 - 
           km_int_T2_HT2*AintL_T2 - km_int_T2_T2_Phase_II*AintL_T2 - km_int_T2_T2_Phase_I*AintL_T2          # Intestinal lumen
         dAint_metab_T2 <- km_sto_T2_HT2*AstoL_T2 + km_sto_T2_T2_Phase_II*AstoL_T2 + km_sto_T2_T2_Phase_I*AstoL_T2 +
           km_int_T2_HT2*AintL_T2 + km_int_T2_T2_Phase_II*AintL_T2 + km_int_T2_T2_Phase_I*AintL_T2  # Amount metabolized in stomach and intestinal lumen
         
         dAfec_T2 <- kfec*AintL_T2	                                 # Feces
         
         Fabs_T2 <- Afec_T2/(Aoral_total_T2+1e-10)
         
         Csto_T2 <- Asto_T2/Vsto
         dAsto_T2 <- ka_sto_T2*AstoL_T2 + Qsto*(Ca_T2 - Csto_T2/Psto_T2)   # Stomach
         
         Cint_T2 <- Aint_T2/Vint
         dAint_T2 <- ka_int_T2*AintL_T2 + Qint*(Ca_T2 - Cint_T2/Pint_T2)        # Intestine
         
         ## IV: see venous compartment below
         Aiv_T2 <- 0 # X.iv_T2[t/tstep + 1]
         
         
         ## intake 
         intake_T2 <- Aair_T2 + Aoral_T2 + As_T2 + Aiv_T2                # Record intake to keep track of mass balance
         
         
         # Storage compartments
         
         Chrt_T2 <- Ahrt_T2/Vhrt
         dAhrt_T2 <- Qhrt*(Ca_T2 - Chrt_T2/Phrt_T2)               # Heart
         Cbrn_T2 <- Abrn_T2/Vbrn
         dAbrn_T2 <- Qbrn*(Ca_T2 - Cbrn_T2/Pbrn_T2)               # Brain
         Cthy_T2 <- Athy_T2/Vthy
         dAthy_T2 <- Qthy*(Ca_T2 - Cthy_T2/Pthy_T2)               # Thymus 
         Cadp_T2 <- Aadp_T2/Vadp
         dAadp_T2 <- Qadp*(Ca_T2 - Cadp_T2/Padp_T2)               # Adipose tissue + yellow marrow
         Cmus_T2 <- Amus_T2/Vmus
         dAmus_T2 <- Qmus*(Ca_T2 - Cmus_T2/Pmus_T2)               # Skeletal muscle
         Cbon_T2 <- Abon_T2/Vbon
         dAbon_T2 <- Qbon*(Ca_T2 - Cbon_T2/Pbon_T2)               # Bone + red marrow
         Cgon_T2 <- Agon_T2/Vgon
         dAgon_T2 <- Qgon*(Ca_T2 - Cgon_T2/Pgon_T2)               # Gonads
         Cpan_T2 <- Apan_T2/Vpan
         dApan_T2 <- Qpan*(Ca_T2 - Cpan_T2/Ppan_T2)               # Pancreas
         Cspl_T2 <- Aspl_T2/Vspl
         dAspl_T2 <- Qspl*(Ca_T2 - Cspl_T2/Pspl_T2)               # Spleen
         Crem_T2 <- Arem_T2/Vrem
         dArem_T2 <- Qrem*(Ca_T2 - Crem_T2/Prem_T2)               # Remaining tissues
         
         # Clearing compartments
         
         # Liver
         
         Qport <- Qsto + Qint + Qpan + Qspl                      # Portal blood flow
         Cport_T2 <- (Qsto*(Csto_T2/Psto_T2) +  # Portal blood 
                     Qint*(Cint_T2/Pint_T2) +                            # concentration
                     Qpan*(Cpan_T2/Ppan_T2) + 
                     Qspl*(Cspl_T2/Pspl_T2))/Qport  
         Cliv_T2 <- Aliv_T2/Vliv  
         dAliv_T2 <- Qliv*Ca_T2 + Qport*Cport_T2 -                       # Liver 
           (Qliv + Qport)*(Cliv_T2/Pliv_T2) -      
           CLliv_T2_HT2*(Cliv_T2/Pliv_T2) -
           CLliv_T2_T2_Phase_II*(Cliv_T2/Pliv_T2) -
           CLliv_T2_T2_Phase_I*(Cliv_T2/Pliv_T2)
         
         dAlivCL_T2 <- CLliv_T2_HT2*(Cliv_T2/Pliv_T2) +
           CLliv_T2_T2_Phase_II*(Cliv_T2/Pliv_T2) +
           CLliv_T2_T2_Phase_I*(Cliv_T2/Pliv_T2)
         # or: Vmaxliv*(Cliv/Pliv)/(kMliv + Cliv/Pliv)     
         
         # Kidneys
         
         Ckid_T2 <- Akid_T2/Vkid
         dAkid_T2 <- Qkid*(Ca_T2 - Ckid_T2/Pkid_T2) - CLkid_T2*Ca_T2        # Kidneys
         dAkidCL_T2 <- CLkid_T2*Ca_T2
         
         # Arterial and venous blood
         
         # Note: 5% of Qc flows through non-included tissues (lymph, a.o.) 
         # TODO: check balance Q
         
         dAa_T2 <- ifelse(include_lung==TRUE, 
                       Qc*(Clun_T2/Plunb_T2),
                       Qc*(Calv_T2/Pab_T2)) -                               # Arterial blood
           (QskE + QskU + Qhrt + Qbrn + Qthy + Qadp +
              Qmus + Qbon + Qgon + Qrem + Qsto + Qint + 
              Qpan + Qspl + Qliv + Qkid + 0.05*Qc)*Ca_T2
         
         dAv_T2 <- Aiv_T2 +                                           # Venous blood
           QskE*CskE_T2/Pskn_T2 +     
           QskU*CskU_T2/Pskn_T2 + 
           Qhrt*Chrt_T2/Phrt_T2 + 
           Qbrn*Cbrn_T2/Pbrn_T2 + 
           Qthy*Cthy_T2/Pthy_T2 + 
           Qadp*Cadp_T2/Padp_T2 + 
           Qmus*Cmus_T2/Pmus_T2 + 
           Qbon*Cbon_T2/Pbon_T2 + 
           Qgon*Cgon_T2/Pgon_T2 + 
           Qrem*Crem_T2/Prem_T2 + 
           (Qliv + Qport)*(Cliv_T2/Pliv_T2) +
           Qkid*Ckid_T2/Pkid_T2 +
           0.05*Qc*Ca_T2 - Qc*Cv_T2                    
         
         # Mass balance
         dCumulativeIntake_T2 <- intake_T2
         MassBal_T2 <- (Aalv_T2 + Alun_T2 + Ass_T2 + AskE_T2 + AskU_T2 + 
                          AstoL_T2 + AintL_T2 + Aint_metab_T2 + Afec_T2 + Asto_T2 + Aint_T2 + 
                          Ahrt_T2 + Abrn_T2 + Athy_T2 + Aadp_T2 + Amus_T2 + Abon_T2 + Agon_T2 + Apan_T2 + Aspl_T2 + Arem_T2 + 
                          Aliv_T2 + AlivCL_T2 + Akid_T2 + AkidCL_T2 + Aa_T2 + Av_T2)#/CumulativeIntake_T2 
         
         #### HT2 ----
         Cair_HT2 <- 0 # X.air_HT2[t/tstep + 1]
         Calv_HT2 <- Aalv_HT2/Valv
         
         Cv_HT2 <- Av_HT2/Vv
         Ca_HT2 <- Aa_HT2/Va
         
         if (include_lung == TRUE) {
           Clun_HT2 <- Alun_HT2/Vlun
           dAalv_HT2 <- Qp*(Cair_HT2 - Calv_HT2) - PAalv_HT2*(Calv_HT2 - Clun_HT2/Pluna_HT2)     # Alveolar air   
           dAlun_HT2 <- Qc*(Cv_HT2 - Clun_HT2/Plunb_HT2) + PAalv_HT2*(Calv_HT2 - Clun_HT2/Pluna_HT2) # Lungs
           Aair_HT2 <- Qp*(Cair_HT2 - Calv_HT2) # outside air
         }
         
         else {
           Clun_HT2 <- Aalv_HT2/Valv
           dAalv_HT2 <- Qp*(Cair_HT2 - Calv_HT2) + Qc*(Cv_HT2 - Calv_HT2/Pab_HT2)  # Alveolar air
           dAlun_HT2 <- 0 # excluded
           Aair_HT2 <- Qp*(Cair_HT2 - Calv_HT2) # outside air
         }
         
         ## Skin contact
         Css_HT2 <- Ass_HT2/Vss # Vss in L (dm3)
         
         As_HT2 <- 0 # X.skin_HT2[t/tstep + 1] # amount of substance added to the skin at time t
         dAss_HT2  <- As_HT2 - Ks_HT2*fss*BSA*Css_HT2/100                             # Skin surface; divided by 100 to correct for BSA unit (m2 -> dm2)
         
         CskE_HT2 <- ifelse(fss!=0, AskE_HT2/(fss*Vskn), 0)
         dAskE_HT2 <- Ks_HT2*fss*BSA*Css_HT2/100 + QskE*(Ca_HT2 - CskE_HT2/Pskn_HT2)      # Exposed skin
         
         CskU_HT2 <- AskU_HT2/((1 - fss)*Vskn)
         dAskU_HT2 <- QskU*(Ca_HT2 - CskU_HT2/Pskn_HT2)                       # Unexposed skin
         
         ## Oral
         #Aoral_HT2 <- X.oral_HT2[t/tstep + 1] 
         Aoral_HT2 <- X.oral_HT2(t) 
         
         dAoral_total_HT2 <- Aoral_HT2
         dAstoL_HT2 <- Aoral_HT2 - ka_sto_HT2*AstoL_HT2 - kse*AstoL_HT2 + km_sto_T2_HT2*AstoL_T2 - 
           km_sto_HT2_HT2_Phase_II*AstoL_HT2 - km_sto_HT2_HT2_Phase_I*AstoL_HT2                # Stomach lumen
         dAintL_HT2 <- kse*AstoL_HT2 - ka_int_HT2*AintL_HT2 - kfec*AintL_HT2 + km_int_T2_HT2*AintL_T2 - 
           km_int_HT2_HT2_Phase_II*AintL_HT2 - km_int_HT2_HT2_Phase_I*AintL_HT2          # Intestinal lumen
         dAint_metab_HT2 <- km_sto_HT2_HT2_Phase_II*AstoL_HT2 + km_sto_HT2_HT2_Phase_I*AstoL_HT2 +
           km_int_HT2_HT2_Phase_II*AintL_HT2 + km_int_HT2_HT2_Phase_I*AintL_HT2 # Amount metabolized in stomach and intestinal lumen
         
         dAfec_HT2 <- kfec*AintL_HT2	                                 # Feces
         
         Fabs_HT2 <- Afec_HT2/(Aoral_total_HT2+1e-10)
         
         Csto_HT2 <- Asto_HT2/Vsto
         dAsto_HT2 <- ka_sto_HT2*AstoL_HT2 + Qsto*(Ca_HT2 - Csto_HT2/Psto_HT2)   # Stomach
         
         Cint_HT2 <- Aint_HT2/Vint
         dAint_HT2 <- ka_int_HT2*AintL_HT2 + Qint*(Ca_HT2 - Cint_HT2/Pint_HT2)        # Intestine
         
         ## IV: see venous compartment below
         Aiv_HT2 <- 0 # X.iv_HT2[t/tstep + 1]
         
         
         ## intake 
         intake_HT2 <- Aair_HT2 + Aoral_HT2 + As_HT2 + Aiv_HT2                # Record intake to keep track of mass balance
         
         
         # Storage compartments
         
         Chrt_HT2 <- Ahrt_HT2/Vhrt
         dAhrt_HT2 <- Qhrt*(Ca_HT2 - Chrt_HT2/Phrt_HT2)               # Heart
         Cbrn_HT2 <- Abrn_HT2/Vbrn
         dAbrn_HT2 <- Qbrn*(Ca_HT2 - Cbrn_HT2/Pbrn_HT2)               # Brain
         Cthy_HT2 <- Athy_HT2/Vthy
         dAthy_HT2 <- Qthy*(Ca_HT2 - Cthy_HT2/Pthy_HT2)               # Thymus 
         Cadp_HT2 <- Aadp_HT2/Vadp
         dAadp_HT2 <- Qadp*(Ca_HT2 - Cadp_HT2/Padp_HT2)               # Adipose tissue + yellow marrow
         Cmus_HT2 <- Amus_HT2/Vmus
         dAmus_HT2 <- Qmus*(Ca_HT2 - Cmus_HT2/Pmus_HT2)               # Skeletal muscle
         Cbon_HT2 <- Abon_HT2/Vbon
         dAbon_HT2 <- Qbon*(Ca_HT2 - Cbon_HT2/Pbon_HT2)               # Bone + red marrow
         Cgon_HT2 <- Agon_HT2/Vgon
         dAgon_HT2 <- Qgon*(Ca_HT2 - Cgon_HT2/Pgon_HT2)               # Gonads
         Cpan_HT2 <- Apan_HT2/Vpan
         dApan_HT2 <- Qpan*(Ca_HT2 - Cpan_HT2/Ppan_HT2)               # Pancreas
         Cspl_HT2 <- Aspl_HT2/Vspl
         dAspl_HT2 <- Qspl*(Ca_HT2 - Cspl_HT2/Pspl_HT2)               # Spleen
         Crem_HT2 <- Arem_HT2/Vrem
         dArem_HT2 <- Qrem*(Ca_HT2 - Crem_HT2/Prem_HT2)               # Remaining tissues
         
         # Clearing compartments
         
         # Liver
         
         Cport_HT2 <- (Qsto*(Csto_HT2/Psto_HT2) + # Portal blood 
                        Qint*(Cint_HT2/Pint_HT2) +                            # concentration
                        Qpan*(Cpan_HT2/Ppan_HT2) + 
                        Qspl*(Cspl_HT2/Pspl_HT2))/Qport  
         Cliv_HT2 <- Aliv_HT2/Vliv  
         dAliv_HT2 <- Qliv*Ca_HT2 + Qport*Cport_HT2 -                       # Liver 
           (Qliv + Qport)*(Cliv_HT2/Pliv_HT2) -      
           CLliv_HT2_HT2_Phase_II*(Cliv_HT2/Pliv_HT2) -
           CLliv_HT2_HT2_Phase_I*(Cliv_HT2/Pliv_HT2) +
           CLliv_T2_HT2*(Cliv_T2/Pliv_T2) # metabolic conversion of T2 into HT2
         
         dAlivCL_HT2 <- CLliv_HT2_HT2_Phase_II*(Cliv_HT2/Pliv_HT2) +
           CLliv_HT2_HT2_Phase_I*(Cliv_HT2/Pliv_HT2)
         # or: Vmaxliv*(Cliv/Pliv)/(kMliv + Cliv/Pliv)     
         
         # Kidneys
         
         Ckid_HT2 <- Akid_HT2/Vkid
         dAkid_HT2 <- Qkid*(Ca_HT2 - Ckid_HT2/Pkid_HT2) - CLkid_HT2*Ca_HT2        # Kidneys
         dAkidCL_HT2 <- CLkid_HT2*Ca_HT2
         
         # Arterial and venous blood
         
         # Note: 5% of Qc flows through non-included tissues (lymph, a.o.) 
         # TODO: check balance Q
         
         dAa_HT2 <- ifelse(include_lung==TRUE, 
                          Qc*(Clun_HT2/Plunb_HT2),
                          Qc*(Calv_HT2/Pab_HT2)) -                               # Arterial blood
           (QskE + QskU + Qhrt + Qbrn + Qthy + Qadp +
              Qmus + Qbon + Qgon + Qrem + Qsto + Qint + 
              Qpan + Qspl + Qliv + Qkid + 0.05*Qc)*Ca_HT2
         
         dAv_HT2 <- Aiv_HT2 +                                           # Venous blood
           QskE*CskE_HT2/Pskn_HT2 + 
           QskU*CskU_HT2/Pskn_HT2 + 
           Qhrt*Chrt_HT2/Phrt_HT2 + 
           Qbrn*Cbrn_HT2/Pbrn_HT2 + 
           Qthy*Cthy_HT2/Pthy_HT2 + 
           Qadp*Cadp_HT2/Padp_HT2 + 
           Qmus*Cmus_HT2/Pmus_HT2 + 
           Qbon*Cbon_HT2/Pbon_HT2 + 
           Qgon*Cgon_HT2/Pgon_HT2 + 
           Qrem*Crem_HT2/Prem_HT2 + 
           (Qliv + Qport)*(Cliv_HT2/Pliv_HT2) +
           Qkid*Ckid_HT2/Pkid_HT2 +
           0.05*Qc*Ca_HT2 - Qc*Cv_HT2                    
         
         # Mass balance
         dCumulativeIntake_HT2 <- intake_HT2 #+ km_sto_T2*AstoL_T2 + km_int_T2*AstoL_T2 + CLliv_T2*(Cliv_T2/Pliv_T2) 
         MassBal_HT2 <- (Aalv_HT2 + Alun_HT2 + Ass_HT2 + AskE_HT2 + AskU_HT2 + AstoL_HT2 + 
                          AintL_HT2 + Aint_metab_HT2 + Asto_HT2 + Aint_HT2 + Afec_HT2 + Ahrt_HT2 + Abrn_HT2 + Athy_HT2 + Aadp_HT2 + Amus_HT2 + Abon_HT2 + 
                          Agon_HT2 + Arem_HT2 + Apan_HT2 + Aspl_HT2 + Aliv_HT2 + AlivCL_HT2 + Akid_HT2 + AkidCL_HT2 + Aa_HT2 + Av_HT2)#/CumulativeIntake_HT2 
         
         #### T2_Phase_I ----
         Cair_T2_Phase_I <- 0 # X.air_T2_Phase_I[t/tstep + 1]
         Calv_T2_Phase_I <- Aalv_T2_Phase_I/Valv
         
         Cv_T2_Phase_I <- Av_T2_Phase_I/Vv
         Ca_T2_Phase_I <- Aa_T2_Phase_I/Va
         
         if (include_lung == TRUE) {
           Clun_T2_Phase_I <- Alun_T2_Phase_I/Vlun
           dAalv_T2_Phase_I <- Qp*(Cair_T2_Phase_I - Calv_T2_Phase_I) - PAalv_T2_Phase_I*(Calv_T2_Phase_I - Clun_T2_Phase_I/Pluna_T2_Phase_I)     # Alveolar air   
           dAlun_T2_Phase_I <- Qc*(Cv_T2_Phase_I - Clun_T2_Phase_I/Plunb_T2_Phase_I) + PAalv_T2_Phase_I*(Calv_T2_Phase_I - Clun_T2_Phase_I/Pluna_T2_Phase_I) # Lungs
           Aair_T2_Phase_I <- Qp*(Cair_T2_Phase_I - Calv_T2_Phase_I) # outside air
         }
         
         else {
           Clun_T2_Phase_I <- Aalv_T2_Phase_I/Valv
           dAalv_T2_Phase_I <- Qp*(Cair_T2_Phase_I - Calv_T2_Phase_I) + Qc*(Cv_T2_Phase_I - Calv_T2_Phase_I/Pab_T2_Phase_I)  # Alveolar air
           dAlun_T2_Phase_I <- 0 # excluded
           Aair_T2_Phase_I <- Qp*(Cair_T2_Phase_I - Calv_T2_Phase_I) # outside air
         }
         
         ## Skin contact
         Css_T2_Phase_I <- Ass_T2_Phase_I/Vss # Vss in L (dm3)
         
         As_T2_Phase_I <- 0 # X.skin_T2_Phase_I[t/tstep + 1] # amount of substance added to the skin at time t
         dAss_T2_Phase_I  <- As_T2_Phase_I - Ks_T2_Phase_I*fss*BSA*Css_T2_Phase_I/100                             # Skin surface; divided by 100 to correct for BSA unit (m2 -> dm2)
         
         CskE_T2_Phase_I <- ifelse(fss!=0, AskE_T2_Phase_I/(fss*Vskn), 0)
         dAskE_T2_Phase_I <- Ks_T2_Phase_I*fss*BSA*Css_T2_Phase_I/100 + QskE*(Ca_T2_Phase_I - CskE_T2_Phase_I/Pskn_T2_Phase_I)      # Exposed skin
         
         CskU_T2_Phase_I <- AskU_T2_Phase_I/((1 - fss)*Vskn)
         dAskU_T2_Phase_I <- QskU*(Ca_T2_Phase_I - CskU_T2_Phase_I/Pskn_T2_Phase_I)                       # Unexposed skin
         
         ## Oral
         #Aoral_T2_Phase_I <- X.oral_T2_Phase_I[t/tstep + 1] 
         Aoral_T2_Phase_I <- 0 # X.oral_T2_Phase_I(t) 
         
         dAoral_total_T2_Phase_I <- Aoral_T2_Phase_I
         dAstoL_T2_Phase_I <- Aoral_T2_Phase_I - ka_sto_T2_Phase_I*AstoL_T2_Phase_I - kse*AstoL_T2_Phase_I + km_sto_T2_T2_Phase_I*AstoL_T2 - 
           km_sto_T2_Phase_I*AstoL_T2_Phase_I                # Stomach lumen
         dAintL_T2_Phase_I <- kse*AstoL_T2_Phase_I - ka_int_T2_Phase_I*AintL_T2_Phase_I - kfec*AintL_T2_Phase_I + km_int_T2_T2_Phase_I*AintL_T2 - 
           km_int_T2_Phase_I*AintL_T2_Phase_I          # Intestinal lumen
         dAint_metab_T2_Phase_I <- km_sto_T2_Phase_I*AstoL_T2_Phase_I + km_int_T2_Phase_I*AintL_T2_Phase_I # Amount metabolized in stomach and intestinal lumen
         
         dAfec_T2_Phase_I <- kfec*AintL_T2_Phase_I	                                 # Feces
         
         Fabs_T2_Phase_I <- Afec_T2_Phase_I/(Aoral_total_T2_Phase_I+1e-10)
         
         Csto_T2_Phase_I <- Asto_T2_Phase_I/Vsto
         dAsto_T2_Phase_I <- ka_sto_T2_Phase_I*AstoL_T2_Phase_I + Qsto*(Ca_T2_Phase_I - Csto_T2_Phase_I/Psto_T2_Phase_I)   # Stomach
         
         Cint_T2_Phase_I <- Aint_T2_Phase_I/Vint
         dAint_T2_Phase_I <- ka_int_T2_Phase_I*AintL_T2_Phase_I + Qint*(Ca_T2_Phase_I - Cint_T2_Phase_I/Pint_T2_Phase_I)        # Intestine
         
         ## IV: see venous compartment below
         Aiv_T2_Phase_I <- 0 # X.iv_T2_Phase_I[t/tstep + 1]
         
         
         ## intake 
         intake_T2_Phase_I <- Aair_T2_Phase_I + Aoral_T2_Phase_I + As_T2_Phase_I + Aiv_T2_Phase_I                # Record intake to keep track of mass balance
         
         
         # Storage compartments
         
         Chrt_T2_Phase_I <- Ahrt_T2_Phase_I/Vhrt
         dAhrt_T2_Phase_I <- Qhrt*(Ca_T2_Phase_I - Chrt_T2_Phase_I/Phrt_T2_Phase_I)               # Heart
         Cbrn_T2_Phase_I <- Abrn_T2_Phase_I/Vbrn
         dAbrn_T2_Phase_I <- Qbrn*(Ca_T2_Phase_I - Cbrn_T2_Phase_I/Pbrn_T2_Phase_I)               # Brain
         Cthy_T2_Phase_I <- Athy_T2_Phase_I/Vthy
         dAthy_T2_Phase_I <- Qthy*(Ca_T2_Phase_I - Cthy_T2_Phase_I/Pthy_T2_Phase_I)               # Thymus 
         Cadp_T2_Phase_I <- Aadp_T2_Phase_I/Vadp
         dAadp_T2_Phase_I <- Qadp*(Ca_T2_Phase_I - Cadp_T2_Phase_I/Padp_T2_Phase_I)               # Adipose tissue + yellow marrow
         Cmus_T2_Phase_I <- Amus_T2_Phase_I/Vmus
         dAmus_T2_Phase_I <- Qmus*(Ca_T2_Phase_I - Cmus_T2_Phase_I/Pmus_T2_Phase_I)               # Skeletal muscle
         Cbon_T2_Phase_I <- Abon_T2_Phase_I/Vbon
         dAbon_T2_Phase_I <- Qbon*(Ca_T2_Phase_I - Cbon_T2_Phase_I/Pbon_T2_Phase_I)               # Bone + red marrow
         Cgon_T2_Phase_I <- Agon_T2_Phase_I/Vgon
         dAgon_T2_Phase_I <- Qgon*(Ca_T2_Phase_I - Cgon_T2_Phase_I/Pgon_T2_Phase_I)               # Gonads
         Cpan_T2_Phase_I <- Apan_T2_Phase_I/Vpan
         dApan_T2_Phase_I <- Qpan*(Ca_T2_Phase_I - Cpan_T2_Phase_I/Ppan_T2_Phase_I)               # Pancreas
         Cspl_T2_Phase_I <- Aspl_T2_Phase_I/Vspl
         dAspl_T2_Phase_I <- Qspl*(Ca_T2_Phase_I - Cspl_T2_Phase_I/Pspl_T2_Phase_I)               # Spleen
         Crem_T2_Phase_I <- Arem_T2_Phase_I/Vrem
         dArem_T2_Phase_I <- Qrem*(Ca_T2_Phase_I - Crem_T2_Phase_I/Prem_T2_Phase_I)               # Remaining tissues
         
         # Clearing compartments
         
         # Liver
         
         Cport_T2_Phase_I <- (Qsto*(Csto_T2_Phase_I/Psto_T2_Phase_I) +  # Portal blood 
                             Qint*(Cint_T2_Phase_I/Pint_T2_Phase_I) +                            # concentration
                             Qpan*(Cpan_T2_Phase_I/Ppan_T2_Phase_I) +  
                             Qspl*(Cspl_T2_Phase_I/Pspl_T2_Phase_I))/Qport  
         Cliv_T2_Phase_I <- Aliv_T2_Phase_I/Vliv  
         dAliv_T2_Phase_I <- Qliv*Ca_T2_Phase_I + Qport*Cport_T2_Phase_I -                       # Liver 
           (Qliv + Qport)*(Cliv_T2_Phase_I/Pliv_T2_Phase_I) -      
           CLliv_T2_Phase_I*(Cliv_T2_Phase_I/Pliv_T2_Phase_I) +
           CLliv_T2_T2_Phase_I*(Cliv_T2/Pliv_T2) # metabolic conversion of T2 into T2_Phase_I
         
         dAlivCL_T2_Phase_I <- CLliv_T2_Phase_I*(Cliv_T2_Phase_I/Pliv_T2_Phase_I)
         # or: Vmaxliv*(Cliv/Pliv)/(kMliv + Cliv/Pliv)     
         
         # Kidneys
         
         Ckid_T2_Phase_I <- Akid_T2_Phase_I/Vkid
         dAkid_T2_Phase_I <- Qkid*(Ca_T2_Phase_I - Ckid_T2_Phase_I/Pkid_T2_Phase_I) - CLkid_T2_Phase_I*Ca_T2_Phase_I        # Kidneys
         dAkidCL_T2_Phase_I <- CLkid_T2_Phase_I*Ca_T2_Phase_I
         
         # Arterial and venous blood
         
         # Note: 5% of Qc flows through non-included tissues (lymph, a.o.) 
         # TODO: check balance Q
         
         dAa_T2_Phase_I <- ifelse(include_lung==TRUE, 
                               Qc*(Clun_T2_Phase_I/Plunb_T2_Phase_I),
                               Qc*(Calv_T2_Phase_I/Pab_T2_Phase_I)) -                               # Arterial blood
           (QskE + QskU + Qhrt + Qbrn + Qthy + Qadp +
              Qmus + Qbon + Qgon + Qrem + Qsto + Qint + 
              Qpan + Qspl + Qliv + Qkid + 0.05*Qc)*Ca_T2_Phase_I
         
         dAv_T2_Phase_I <- Aiv_T2_Phase_I +                                           # Venous blood
           QskE*CskE_T2_Phase_I/Pskn_T2_Phase_I + 
           QskU*CskU_T2_Phase_I/Pskn_T2_Phase_I + 
           Qhrt*Chrt_T2_Phase_I/Phrt_T2_Phase_I + 
           Qbrn*Cbrn_T2_Phase_I/Pbrn_T2_Phase_I + 
           Qthy*Cthy_T2_Phase_I/Pthy_T2_Phase_I + 
           Qadp*Cadp_T2_Phase_I/Padp_T2_Phase_I + 
           Qmus*Cmus_T2_Phase_I/Pmus_T2_Phase_I + 
           Qbon*Cbon_T2_Phase_I/Pbon_T2_Phase_I + 
           Qgon*Cgon_T2_Phase_I/Pgon_T2_Phase_I + 
           Qrem*Crem_T2_Phase_I/Prem_T2_Phase_I + 
           (Qliv + Qport)*(Cliv_T2_Phase_I/Pliv_T2_Phase_I) +
           Qkid*Ckid_T2_Phase_I/Pkid_T2_Phase_I +
           0.05*Qc*Ca_T2_Phase_I - Qc*Cv_T2_Phase_I                    
         
         # Mass balance
         dCumulativeIntake_T2_Phase_I <- intake_T2_Phase_I #+ km_sto_T2*AstoL_T2 + km_int_T2*AstoL_T2 + CLliv_T2*(Cliv_T2/Pliv_T2) 
         MassBal_T2_Phase_I <- (Aalv_T2_Phase_I + Alun_T2_Phase_I + Ass_T2_Phase_I + AskE_T2_Phase_I + AskU_T2_Phase_I + AstoL_T2_Phase_I + 
                               AintL_T2_Phase_I + Aint_metab_T2_Phase_I + Asto_T2_Phase_I + Aint_T2_Phase_I + Afec_T2_Phase_I + Ahrt_T2_Phase_I + Abrn_T2_Phase_I + Athy_T2_Phase_I + Aadp_T2_Phase_I + Amus_T2_Phase_I + Abon_T2_Phase_I + 
                               Agon_T2_Phase_I + Arem_T2_Phase_I + Apan_T2_Phase_I + Aspl_T2_Phase_I + Aliv_T2_Phase_I + AlivCL_T2_Phase_I + Akid_T2_Phase_I + AkidCL_T2_Phase_I + Aa_T2_Phase_I + Av_T2_Phase_I)#/CumulativeIntake_T2_Phase_I 
         
         #### HT2_Phase_I ----
         Cair_HT2_Phase_I <- 0 # X.air_HT2_Phase_I[t/tstep + 1]
         Calv_HT2_Phase_I <- Aalv_HT2_Phase_I/Valv
         
         Cv_HT2_Phase_I <- Av_HT2_Phase_I/Vv
         Ca_HT2_Phase_I <- Aa_HT2_Phase_I/Va
         
         if (include_lung == TRUE) {
           Clun_HT2_Phase_I <- Alun_HT2_Phase_I/Vlun
           dAalv_HT2_Phase_I <- Qp*(Cair_HT2_Phase_I - Calv_HT2_Phase_I) - PAalv_HT2_Phase_I*(Calv_HT2_Phase_I - Clun_HT2_Phase_I/Pluna_HT2_Phase_I)     # Alveolar air   
           dAlun_HT2_Phase_I <- Qc*(Cv_HT2_Phase_I - Clun_HT2_Phase_I/Plunb_HT2_Phase_I) + PAalv_HT2_Phase_I*(Calv_HT2_Phase_I - Clun_HT2_Phase_I/Pluna_HT2_Phase_I) # Lungs
           Aair_HT2_Phase_I <- Qp*(Cair_HT2_Phase_I - Calv_HT2_Phase_I) # outside air
         }
         
         else {
           Clun_HT2_Phase_I <- Aalv_HT2_Phase_I/Valv
           dAalv_HT2_Phase_I <- Qp*(Cair_HT2_Phase_I - Calv_HT2_Phase_I) + Qc*(Cv_HT2_Phase_I - Calv_HT2_Phase_I/Pab_HT2_Phase_I)  # Alveolar air
           dAlun_HT2_Phase_I <- 0 # excluded
           Aair_HT2_Phase_I <- Qp*(Cair_HT2_Phase_I - Calv_HT2_Phase_I) # outside air
         }
         
         ## Skin contact
         Css_HT2_Phase_I <- Ass_HT2_Phase_I/Vss # Vss in L (dm3)
         
         As_HT2_Phase_I <- 0 # X.skin_HT2_Phase_I[t/tstep + 1] # amount of substance added to the skin at time t
         dAss_HT2_Phase_I  <- As_HT2_Phase_I - Ks_HT2_Phase_I*fss*BSA*Css_HT2_Phase_I/100                             # Skin surface; divided by 100 to correct for BSA unit (m2 -> dm2)
         
         CskE_HT2_Phase_I <- ifelse(fss!=0, AskE_HT2_Phase_I/(fss*Vskn), 0)
         dAskE_HT2_Phase_I <- Ks_HT2_Phase_I*fss*BSA*Css_HT2_Phase_I/100 + QskE*(Ca_HT2_Phase_I - CskE_HT2_Phase_I/Pskn_HT2_Phase_I)      # Exposed skin
         
         CskU_HT2_Phase_I <- AskU_HT2_Phase_I/((1 - fss)*Vskn)
         dAskU_HT2_Phase_I <- QskU*(Ca_HT2_Phase_I - CskU_HT2_Phase_I/Pskn_HT2_Phase_I)                       # Unexposed skin
         
         ## Oral
         #Aoral_HT2_Phase_I <- X.oral_HT2_Phase_I[t/tstep + 1] 
         Aoral_HT2_Phase_I <- 0 # X.oral_HT2_Phase_I(t) 
         
         dAoral_total_HT2_Phase_I <- Aoral_HT2_Phase_I
         dAstoL_HT2_Phase_I <- Aoral_HT2_Phase_I - ka_sto_HT2_Phase_I*AstoL_HT2_Phase_I - kse*AstoL_HT2_Phase_I + km_sto_HT2_HT2_Phase_I*AstoL_HT2  # Stomach lumen
         dAintL_HT2_Phase_I <- kse*AstoL_HT2_Phase_I - ka_int_HT2_Phase_I*AintL_HT2_Phase_I - kfec*AintL_HT2_Phase_I + km_int_HT2_HT2_Phase_I*AintL_HT2  # Intestinal lumen
         dAint_metab_HT2_Phase_I <- 0 # Amount metabolized in stomach and intestinal lumen
         
         dAfec_HT2_Phase_I <- kfec*AintL_HT2_Phase_I	                                 # Feces
         
         Fabs_HT2_Phase_I <- Afec_HT2_Phase_I/(Aoral_total_HT2_Phase_I+1e-10)
         
         Csto_HT2_Phase_I <- Asto_HT2_Phase_I/Vsto
         dAsto_HT2_Phase_I <- ka_sto_HT2_Phase_I*AstoL_HT2_Phase_I + Qsto*(Ca_HT2_Phase_I - Csto_HT2_Phase_I/Psto_HT2_Phase_I)   # Stomach
         
         Cint_HT2_Phase_I <- Aint_HT2_Phase_I/Vint
         dAint_HT2_Phase_I <- ka_int_HT2_Phase_I*AintL_HT2_Phase_I + Qint*(Ca_HT2_Phase_I - Cint_HT2_Phase_I/Pint_HT2_Phase_I)        # Intestine
         
         ## IV: see venous compartment below
         Aiv_HT2_Phase_I <- 0 # X.iv_HT2_Phase_I[t/tstep + 1]
         
         
         ## intake 
         intake_HT2_Phase_I <- Aair_HT2_Phase_I + Aoral_HT2_Phase_I + As_HT2_Phase_I + Aiv_HT2_Phase_I                # Record intake to keep track of mass balance
         
         
         # Storage compartments
         
         Chrt_HT2_Phase_I <- Ahrt_HT2_Phase_I/Vhrt
         dAhrt_HT2_Phase_I <- Qhrt*(Ca_HT2_Phase_I - Chrt_HT2_Phase_I/Phrt_HT2_Phase_I)               # Heart
         Cbrn_HT2_Phase_I <- Abrn_HT2_Phase_I/Vbrn
         dAbrn_HT2_Phase_I <- Qbrn*(Ca_HT2_Phase_I - Cbrn_HT2_Phase_I/Pbrn_HT2_Phase_I)               # Brain
         Cthy_HT2_Phase_I <- Athy_HT2_Phase_I/Vthy
         dAthy_HT2_Phase_I <- Qthy*(Ca_HT2_Phase_I - Cthy_HT2_Phase_I/Pthy_HT2_Phase_I)               # Thymus 
         Cadp_HT2_Phase_I <- Aadp_HT2_Phase_I/Vadp
         dAadp_HT2_Phase_I <- Qadp*(Ca_HT2_Phase_I - Cadp_HT2_Phase_I/Padp_HT2_Phase_I)               # Adipose tissue + yellow marrow
         Cmus_HT2_Phase_I <- Amus_HT2_Phase_I/Vmus
         dAmus_HT2_Phase_I <- Qmus*(Ca_HT2_Phase_I - Cmus_HT2_Phase_I/Pmus_HT2_Phase_I)               # Skeletal muscle
         Cbon_HT2_Phase_I <- Abon_HT2_Phase_I/Vbon
         dAbon_HT2_Phase_I <- Qbon*(Ca_HT2_Phase_I - Cbon_HT2_Phase_I/Pbon_HT2_Phase_I)               # Bone + red marrow
         Cgon_HT2_Phase_I <- Agon_HT2_Phase_I/Vgon
         dAgon_HT2_Phase_I <- Qgon*(Ca_HT2_Phase_I - Cgon_HT2_Phase_I/Pgon_HT2_Phase_I)               # Gonads
         Cpan_HT2_Phase_I <- Apan_HT2_Phase_I/Vpan
         dApan_HT2_Phase_I <- Qpan*(Ca_HT2_Phase_I - Cpan_HT2_Phase_I/Ppan_HT2_Phase_I)               # Pancreas
         Cspl_HT2_Phase_I <- Aspl_HT2_Phase_I/Vspl
         dAspl_HT2_Phase_I <- Qspl*(Ca_HT2_Phase_I - Cspl_HT2_Phase_I/Pspl_HT2_Phase_I)               # Spleen
         Crem_HT2_Phase_I <- Arem_HT2_Phase_I/Vrem
         dArem_HT2_Phase_I <- Qrem*(Ca_HT2_Phase_I - Crem_HT2_Phase_I/Prem_HT2_Phase_I)               # Remaining tissues
         
         # Clearing compartments
         
         # Liver
         
         Cport_HT2_Phase_I <- (Qsto*(Csto_HT2_Phase_I/Psto_HT2_Phase_I) +  # Portal blood 
                              Qint*(Cint_HT2_Phase_I/Pint_HT2_Phase_I) +                            # concentration
                              Qpan*(Cpan_HT2_Phase_I/Ppan_HT2_Phase_I) + 
                              Qspl*(Cspl_HT2_Phase_I/Pspl_HT2_Phase_I))/Qport  
         Cliv_HT2_Phase_I <- Aliv_HT2_Phase_I/Vliv  
         dAliv_HT2_Phase_I <- Qliv*Ca_HT2_Phase_I + Qport*Cport_HT2_Phase_I -                       # Liver 
           (Qliv + Qport)*(Cliv_HT2_Phase_I/Pliv_HT2_Phase_I) +
           CLliv_HT2_HT2_Phase_I*(Cliv_HT2/Pliv_HT2) # metabolic conversion of T2 into HT2_Phase_I
         
         dAlivCL_HT2_Phase_I <- 0
         # or: Vmaxliv*(Cliv/Pliv)/(kMliv + Cliv/Pliv)     
         
         # Kidneys
         
         Ckid_HT2_Phase_I <- Akid_HT2_Phase_I/Vkid
         dAkid_HT2_Phase_I <- Qkid*(Ca_HT2_Phase_I - Ckid_HT2_Phase_I/Pkid_HT2_Phase_I) - CLkid_HT2_Phase_I*Ca_HT2_Phase_I        # Kidneys
         dAkidCL_HT2_Phase_I <- CLkid_HT2_Phase_I*Ca_HT2_Phase_I
         
         # Arterial and venous blood
         
         # Note: 5% of Qc flows through non-included tissues (lymph, a.o.) 
         # TODO: check balance Q
         
         dAa_HT2_Phase_I <- ifelse(include_lung==TRUE, 
                                Qc*(Clun_HT2_Phase_I/Plunb_HT2_Phase_I),
                                Qc*(Calv_HT2_Phase_I/Pab_HT2_Phase_I)) -                               # Arterial blood
           (QskE + QskU + Qhrt + Qbrn + Qthy + Qadp +
              Qmus + Qbon + Qgon + Qrem + Qsto + Qint + 
              Qpan + Qspl + Qliv + Qkid + 0.05*Qc)*Ca_HT2_Phase_I
         
         dAv_HT2_Phase_I <- Aiv_HT2_Phase_I +                                           # Venous blood
           QskE*CskE_HT2_Phase_I/Pskn_HT2_Phase_I + 
           QskU*CskU_HT2_Phase_I/Pskn_HT2_Phase_I + 
           Qhrt*Chrt_HT2_Phase_I/Phrt_HT2_Phase_I + 
           Qbrn*Cbrn_HT2_Phase_I/Pbrn_HT2_Phase_I + 
           Qthy*Cthy_HT2_Phase_I/Pthy_HT2_Phase_I + 
           Qadp*Cadp_HT2_Phase_I/Padp_HT2_Phase_I + 
           Qmus*Cmus_HT2_Phase_I/Pmus_HT2_Phase_I + 
           Qbon*Cbon_HT2_Phase_I/Pbon_HT2_Phase_I + 
           Qgon*Cgon_HT2_Phase_I/Pgon_HT2_Phase_I + 
           Qrem*Crem_HT2_Phase_I/Prem_HT2_Phase_I + 
           (Qliv + Qport)*(Cliv_HT2_Phase_I/Pliv_HT2_Phase_I) +
           Qkid*Ckid_HT2_Phase_I/Pkid_HT2_Phase_I +
           0.05*Qc*Ca_HT2_Phase_I - Qc*Cv_HT2_Phase_I                    
         
         # Mass balance
         dCumulativeIntake_HT2_Phase_I <- intake_HT2_Phase_I #+ km_sto_T2*AstoL_T2 + km_int_T2*AstoL_T2 + CLliv_T2*(Cliv_T2/Pliv_T2) 
         MassBal_HT2_Phase_I <- (Aalv_HT2_Phase_I + Alun_HT2_Phase_I + Ass_HT2_Phase_I + AskE_HT2_Phase_I + AskU_HT2_Phase_I + AstoL_HT2_Phase_I + 
                                AintL_HT2_Phase_I + Aint_metab_HT2_Phase_I + Asto_HT2_Phase_I + Aint_HT2_Phase_I + Afec_HT2_Phase_I + Ahrt_HT2_Phase_I + Abrn_HT2_Phase_I + Athy_HT2_Phase_I + Aadp_HT2_Phase_I + Amus_HT2_Phase_I + Abon_HT2_Phase_I + 
                                Agon_HT2_Phase_I + Arem_HT2_Phase_I + Apan_HT2_Phase_I + Aspl_HT2_Phase_I + Aliv_HT2_Phase_I + AlivCL_HT2_Phase_I + Akid_HT2_Phase_I + AkidCL_HT2_Phase_I + Aa_HT2_Phase_I + Av_HT2_Phase_I)#/CumulativeIntake_HT2_Phase_I 
         
         
         #### T2_Phase_II ----
         Cair_T2_Phase_II <- 0 # X.air_T2_Phase_II[t/tstep + 1]
         Calv_T2_Phase_II <- Aalv_T2_Phase_II/Valv
         
         Cv_T2_Phase_II <- Av_T2_Phase_II/Vv
         Ca_T2_Phase_II <- Aa_T2_Phase_II/Va
         
         if (include_lung == TRUE) {
           Clun_T2_Phase_II <- Alun_T2_Phase_II/Vlun
           dAalv_T2_Phase_II <- Qp*(Cair_T2_Phase_II - Calv_T2_Phase_II) - PAalv_T2_Phase_II*(Calv_T2_Phase_II - Clun_T2_Phase_II/Pluna_T2_Phase_II)     # Alveolar air   
           dAlun_T2_Phase_II <- Qc*(Cv_T2_Phase_II - Clun_T2_Phase_II/Plunb_T2_Phase_II) + PAalv_T2_Phase_II*(Calv_T2_Phase_II - Clun_T2_Phase_II/Pluna_T2_Phase_II) # Lungs
           Aair_T2_Phase_II <- Qp*(Cair_T2_Phase_II - Calv_T2_Phase_II) # outside air
         }
         
         else {
           Clun_T2_Phase_II <- Aalv_T2_Phase_II/Valv
           dAalv_T2_Phase_II <- Qp*(Cair_T2_Phase_II - Calv_T2_Phase_II) + Qc*(Cv_T2_Phase_II - Calv_T2_Phase_II/Pab_T2_Phase_II)  # Alveolar air
           dAlun_T2_Phase_II <- 0 # excluded
           Aair_T2_Phase_II <- Qp*(Cair_T2_Phase_II - Calv_T2_Phase_II) # outside air
         }
         
         ## Skin contact
         Css_T2_Phase_II <- Ass_T2_Phase_II/Vss # Vss in L (dm3)
         
         As_T2_Phase_II <- 0 # X.skin_T2_Phase_II[t/tstep + 1] # amount of substance added to the skin at time t
         dAss_T2_Phase_II  <- As_T2_Phase_II - Ks_T2_Phase_II*fss*BSA*Css_T2_Phase_II/100                             # Skin surface; divided by 100 to correct for BSA unit (m2 -> dm2)
         
         CskE_T2_Phase_II <- ifelse(fss!=0, AskE_T2_Phase_II/(fss*Vskn), 0)
         dAskE_T2_Phase_II <- Ks_T2_Phase_II*fss*BSA*Css_T2_Phase_II/100 + QskE*(Ca_T2_Phase_II - CskE_T2_Phase_II/Pskn_T2_Phase_II)      # Exposed skin
         
         CskU_T2_Phase_II <- AskU_T2_Phase_II/((1 - fss)*Vskn)
         dAskU_T2_Phase_II <- QskU*(Ca_T2_Phase_II - CskU_T2_Phase_II/Pskn_T2_Phase_II)                       # Unexposed skin
         
         ## Oral
         #Aoral_T2_Phase_II <- X.oral_T2_Phase_II[t/tstep + 1] 
         Aoral_T2_Phase_II <- 0 # X.oral_T2_Phase_II(t) 
         
         dAoral_total_T2_Phase_II <- Aoral_T2_Phase_II
         dAstoL_T2_Phase_II <- Aoral_T2_Phase_II - ka_sto_T2_Phase_II*AstoL_T2_Phase_II - kse*AstoL_T2_Phase_II + km_sto_T2_T2_Phase_II*AstoL_T2 - 
           km_sto_T2_Phase_II_HT2_Phase_II*AstoL_T2_Phase_II                # Stomach lumen
         dAintL_T2_Phase_II <- kse*AstoL_T2_Phase_II - ka_int_T2_Phase_II*AintL_T2_Phase_II - kfec*AintL_T2_Phase_II + km_int_T2_T2_Phase_II*AintL_T2 - 
           km_int_T2_Phase_II_HT2_Phase_II*AintL_T2_Phase_II          # Intestinal lumen
         dAint_metab_T2_Phase_II <- km_sto_T2_Phase_II_HT2_Phase_II*AstoL_T2_Phase_II +
           km_int_T2_Phase_II_HT2_Phase_II*AintL_T2_Phase_II # Amount metabolized in stomach and intestinal lumen
         
         dAfec_T2_Phase_II <- kfec*AintL_T2_Phase_II	                                 # Feces
         
         Fabs_T2_Phase_II <- Afec_T2_Phase_II/(Aoral_total_T2_Phase_II+1e-10)
         
         Csto_T2_Phase_II <- Asto_T2_Phase_II/Vsto
         dAsto_T2_Phase_II <- ka_sto_T2_Phase_II*AstoL_T2_Phase_II + Qsto*(Ca_T2_Phase_II - Csto_T2_Phase_II/Psto_T2_Phase_II)   # Stomach
         
         Cint_T2_Phase_II <- Aint_T2_Phase_II/Vint
         dAint_T2_Phase_II <- ka_int_T2_Phase_II*AintL_T2_Phase_II + Qint*(Ca_T2_Phase_II - Cint_T2_Phase_II/Pint_T2_Phase_II)        # Intestine
         
         ## IV: see venous compartment below
         Aiv_T2_Phase_II <- 0 # X.iv_T2_Phase_II[t/tstep + 1]
         
         
         ## intake 
         intake_T2_Phase_II <- Aair_T2_Phase_II + Aoral_T2_Phase_II + As_T2_Phase_II + Aiv_T2_Phase_II                # Record intake to keep track of mass balance
         
         
         # Storage compartments
         
         Chrt_T2_Phase_II <- Ahrt_T2_Phase_II/Vhrt
         dAhrt_T2_Phase_II <- Qhrt*(Ca_T2_Phase_II - Chrt_T2_Phase_II/Phrt_T2_Phase_II)               # Heart
         Cbrn_T2_Phase_II <- Abrn_T2_Phase_II/Vbrn
         dAbrn_T2_Phase_II <- Qbrn*(Ca_T2_Phase_II - Cbrn_T2_Phase_II/Pbrn_T2_Phase_II)               # Brain
         Cthy_T2_Phase_II <- Athy_T2_Phase_II/Vthy
         dAthy_T2_Phase_II <- Qthy*(Ca_T2_Phase_II - Cthy_T2_Phase_II/Pthy_T2_Phase_II)               # Thymus 
         Cadp_T2_Phase_II <- Aadp_T2_Phase_II/Vadp
         dAadp_T2_Phase_II <- Qadp*(Ca_T2_Phase_II - Cadp_T2_Phase_II/Padp_T2_Phase_II)               # Adipose tissue + yellow marrow
         Cmus_T2_Phase_II <- Amus_T2_Phase_II/Vmus
         dAmus_T2_Phase_II <- Qmus*(Ca_T2_Phase_II - Cmus_T2_Phase_II/Pmus_T2_Phase_II)               # Skeletal muscle
         Cbon_T2_Phase_II <- Abon_T2_Phase_II/Vbon
         dAbon_T2_Phase_II <- Qbon*(Ca_T2_Phase_II - Cbon_T2_Phase_II/Pbon_T2_Phase_II)               # Bone + red marrow
         Cgon_T2_Phase_II <- Agon_T2_Phase_II/Vgon
         dAgon_T2_Phase_II <- Qgon*(Ca_T2_Phase_II - Cgon_T2_Phase_II/Pgon_T2_Phase_II)               # Gonads
         Cpan_T2_Phase_II <- Apan_T2_Phase_II/Vpan
         dApan_T2_Phase_II <- Qpan*(Ca_T2_Phase_II - Cpan_T2_Phase_II/Ppan_T2_Phase_II)               # Pancreas
         Cspl_T2_Phase_II <- Aspl_T2_Phase_II/Vspl
         dAspl_T2_Phase_II <- Qspl*(Ca_T2_Phase_II - Cspl_T2_Phase_II/Pspl_T2_Phase_II)               # Spleen
         Crem_T2_Phase_II <- Arem_T2_Phase_II/Vrem
         dArem_T2_Phase_II <- Qrem*(Ca_T2_Phase_II - Crem_T2_Phase_II/Prem_T2_Phase_II)               # Remaining tissues
         
         # Clearing compartments
         
         # Liver
         
         Cport_T2_Phase_II <- (Qsto*(Csto_T2_Phase_II/Psto_T2_Phase_II) +  # Portal blood 
                                 Qint*(Cint_T2_Phase_II/Pint_T2_Phase_II) +                            # concentration
                                 Qpan*(Cpan_T2_Phase_II/Ppan_T2_Phase_II) + 
                                 Qspl*(Cspl_T2_Phase_II/Pspl_T2_Phase_II))/Qport  
         
         Cliv_T2_Phase_II <- Aliv_T2_Phase_II/Vliv  
         dAliv_T2_Phase_II <- Qliv*Ca_T2_Phase_II + Qport*Cport_T2_Phase_II -                       # Liver 
           (Qliv + Qport)*(Cliv_T2_Phase_II/Pliv_T2_Phase_II) -      
           CLliv_T2_Phase_II_HT2_Phase_II*(Cliv_T2_Phase_II/Pliv_T2_Phase_II) +
           CLliv_T2_T2_Phase_II*(Cliv_T2/Pliv_T2) # metabolic conversion of T2 into T2_Phase_II
         
         dAlivCL_T2_Phase_II <- CLliv_T2_Phase_II_HT2_Phase_II*(Cliv_T2_Phase_II/Pliv_T2_Phase_II)
         # or: Vmaxliv*(Cliv/Pliv)/(kMliv + Cliv/Pliv)     
         
         # Kidneys
         
         Ckid_T2_Phase_II <- Akid_T2_Phase_II/Vkid
         dAkid_T2_Phase_II <- Qkid*(Ca_T2_Phase_II - Ckid_T2_Phase_II/Pkid_T2_Phase_II) - CLkid_T2_Phase_II*Ca_T2_Phase_II        # Kidneys
         dAkidCL_T2_Phase_II <- CLkid_T2_Phase_II*Ca_T2_Phase_II
         
         # Arterial and venous blood
         
         # Note: 5% of Qc flows through non-included tissues (lymph, a.o.) 
         # TODO: check balance Q
         
         dAa_T2_Phase_II <- ifelse(include_lung==TRUE, 
                                   Qc*(Clun_T2_Phase_II/Plunb_T2_Phase_II),
                                   Qc*(Calv_T2_Phase_II/Pab_T2_Phase_II)) -                               # Arterial blood
           (QskE + QskU + Qhrt + Qbrn + Qthy + Qadp +
              Qmus + Qbon + Qgon + Qrem + Qsto + Qint + 
              Qpan + Qspl + Qliv + Qkid + 0.05*Qc)*Ca_T2_Phase_II
         
         dAv_T2_Phase_II <- Aiv_T2_Phase_II +                                           # Venous blood
           QskE*CskE_T2_Phase_II/Pskn_T2_Phase_II + 
           QskU*CskU_T2_Phase_II/Pskn_T2_Phase_II + 
           Qhrt*Chrt_T2_Phase_II/Phrt_T2_Phase_II + 
           Qbrn*Cbrn_T2_Phase_II/Pbrn_T2_Phase_II + 
           Qthy*Cthy_T2_Phase_II/Pthy_T2_Phase_II + 
           Qadp*Cadp_T2_Phase_II/Padp_T2_Phase_II + 
           Qmus*Cmus_T2_Phase_II/Pmus_T2_Phase_II + 
           Qbon*Cbon_T2_Phase_II/Pbon_T2_Phase_II + 
           Qgon*Cgon_T2_Phase_II/Pgon_T2_Phase_II + 
           Qrem*Crem_T2_Phase_II/Prem_T2_Phase_II + 
           (Qliv + Qport)*(Cliv_T2_Phase_II/Pliv_T2_Phase_II) +
           Qkid*Ckid_T2_Phase_II/Pkid_T2_Phase_II +
           0.05*Qc*Ca_T2_Phase_II - Qc*Cv_T2_Phase_II                    
         
         # Mass balance
         dCumulativeIntake_T2_Phase_II <- intake_T2_Phase_II #+ km_sto_T2*AstoL_T2 + km_int_T2*AstoL_T2 + CLliv_T2*(Cliv_T2/Pliv_T2) 
         MassBal_T2_Phase_II <- (Aalv_T2_Phase_II + Alun_T2_Phase_II + Ass_T2_Phase_II + AskE_T2_Phase_II + AskU_T2_Phase_II + AstoL_T2_Phase_II + 
                                   AintL_T2_Phase_II + Aint_metab_T2_Phase_II + Asto_T2_Phase_II + Aint_T2_Phase_II + Afec_T2_Phase_II + Ahrt_T2_Phase_II + Abrn_T2_Phase_II + Athy_T2_Phase_II + Aadp_T2_Phase_II + Amus_T2_Phase_II + Abon_T2_Phase_II + 
                                   Agon_T2_Phase_II + Arem_T2_Phase_II + Apan_T2_Phase_II + Aspl_T2_Phase_II + Aliv_T2_Phase_II + AlivCL_T2_Phase_II + Akid_T2_Phase_II + AkidCL_T2_Phase_II + Aa_T2_Phase_II + Av_T2_Phase_II)#/CumulativeIntake_T2_Phase_II 
         
         #### HT2_Phase_II ----
         Cair_HT2_Phase_II <- 0 # X.air_HT2_Phase_II[t/tstep + 1]
         Calv_HT2_Phase_II <- Aalv_HT2_Phase_II/Valv
         
         Cv_HT2_Phase_II <- Av_HT2_Phase_II/Vv
         Ca_HT2_Phase_II <- Aa_HT2_Phase_II/Va
         
         if (include_lung == TRUE) {
           Clun_HT2_Phase_II <- Alun_HT2_Phase_II/Vlun
           dAalv_HT2_Phase_II <- Qp*(Cair_HT2_Phase_II - Calv_HT2_Phase_II) - PAalv_HT2_Phase_II*(Calv_HT2_Phase_II - Clun_HT2_Phase_II/Pluna_HT2_Phase_II)     # Alveolar air   
           dAlun_HT2_Phase_II <- Qc*(Cv_HT2_Phase_II - Clun_HT2_Phase_II/Plunb_HT2_Phase_II) + PAalv_HT2_Phase_II*(Calv_HT2_Phase_II - Clun_HT2_Phase_II/Pluna_HT2_Phase_II) # Lungs
           Aair_HT2_Phase_II <- Qp*(Cair_HT2_Phase_II - Calv_HT2_Phase_II) # outside air
         }
         
         else {
           Clun_HT2_Phase_II <- Aalv_HT2_Phase_II/Valv
           dAalv_HT2_Phase_II <- Qp*(Cair_HT2_Phase_II - Calv_HT2_Phase_II) + Qc*(Cv_HT2_Phase_II - Calv_HT2_Phase_II/Pab_HT2_Phase_II)  # Alveolar air
           dAlun_HT2_Phase_II <- 0 # excluded
           Aair_HT2_Phase_II <- Qp*(Cair_HT2_Phase_II - Calv_HT2_Phase_II) # outside air
         }
         
         ## Skin contact
         Css_HT2_Phase_II <- Ass_HT2_Phase_II/Vss # Vss in L (dm3)
         
         As_HT2_Phase_II <- 0 # X.skin_HT2_Phase_II[t/tstep + 1] # amount of substance added to the skin at time t
         dAss_HT2_Phase_II  <- As_HT2_Phase_II - Ks_HT2_Phase_II*fss*BSA*Css_HT2_Phase_II/100                             # Skin surface; divided by 100 to correct for BSA unit (m2 -> dm2)
         
         CskE_HT2_Phase_II <- ifelse(fss!=0, AskE_HT2_Phase_II/(fss*Vskn), 0)
         dAskE_HT2_Phase_II <- Ks_HT2_Phase_II*fss*BSA*Css_HT2_Phase_II/100 + QskE*(Ca_HT2_Phase_II - CskE_HT2_Phase_II/Pskn_HT2_Phase_II)      # Exposed skin
         
         CskU_HT2_Phase_II <- AskU_HT2_Phase_II/((1 - fss)*Vskn)
         dAskU_HT2_Phase_II <- QskU*(Ca_HT2_Phase_II - CskU_HT2_Phase_II/Pskn_HT2_Phase_II)                       # Unexposed skin
         
         ## Oral
         #Aoral_HT2_Phase_II <- X.oral_HT2_Phase_II[t/tstep + 1] 
         Aoral_HT2_Phase_II <- 0 # X.oral_HT2_Phase_II(t) 
         
         dAoral_total_HT2_Phase_II <- Aoral_HT2_Phase_II
         dAstoL_HT2_Phase_II <- Aoral_HT2_Phase_II - ka_sto_HT2_Phase_II*AstoL_HT2_Phase_II - kse*AstoL_HT2_Phase_II + km_sto_HT2_HT2_Phase_II*AstoL_HT2 + km_sto_T2_Phase_II_HT2_Phase_II*AstoL_T2_Phase_II          # Stomach lumen
         dAintL_HT2_Phase_II <- kse*AstoL_HT2_Phase_II - ka_int_HT2_Phase_II*AintL_HT2_Phase_II - kfec*AintL_HT2_Phase_II + km_int_HT2_HT2_Phase_II*AintL_T2 + km_int_T2_Phase_II_HT2_Phase_II*AintL_T2_Phase_II    # Intestinal lumen
         dAint_metab_HT2_Phase_II <- 0 # Amount metabolized in stomach and intestinal lumen
         
         dAfec_HT2_Phase_II <- kfec*AintL_HT2_Phase_II	                                 # Feces
         
         Fabs_HT2_Phase_II <- Afec_HT2_Phase_II/(Aoral_total_HT2_Phase_II+1e-10)
         
         Csto_HT2_Phase_II <- Asto_HT2_Phase_II/Vsto
         dAsto_HT2_Phase_II <- ka_sto_HT2_Phase_II*AstoL_HT2_Phase_II + Qsto*(Ca_HT2_Phase_II - Csto_HT2_Phase_II/Psto_HT2_Phase_II)   # Stomach
         
         Cint_HT2_Phase_II <- Aint_HT2_Phase_II/Vint
         dAint_HT2_Phase_II <- ka_int_HT2_Phase_II*AintL_HT2_Phase_II + Qint*(Ca_HT2_Phase_II - Cint_HT2_Phase_II/Pint_HT2_Phase_II)        # Intestine
         
         ## IV: see venous compartment below
         Aiv_HT2_Phase_II <- 0 # X.iv_HT2_Phase_II[t/tstep + 1]
         
         
         ## intake 
         intake_HT2_Phase_II <- Aair_HT2_Phase_II + Aoral_HT2_Phase_II + As_HT2_Phase_II + Aiv_HT2_Phase_II                # Record intake to keep track of mass balance
         
         
         # Storage compartments
         
         Chrt_HT2_Phase_II <- Ahrt_HT2_Phase_II/Vhrt
         dAhrt_HT2_Phase_II <- Qhrt*(Ca_HT2_Phase_II - Chrt_HT2_Phase_II/Phrt_HT2_Phase_II)               # Heart
         Cbrn_HT2_Phase_II <- Abrn_HT2_Phase_II/Vbrn
         dAbrn_HT2_Phase_II <- Qbrn*(Ca_HT2_Phase_II - Cbrn_HT2_Phase_II/Pbrn_HT2_Phase_II)               # Brain
         Cthy_HT2_Phase_II <- Athy_HT2_Phase_II/Vthy
         dAthy_HT2_Phase_II <- Qthy*(Ca_HT2_Phase_II - Cthy_HT2_Phase_II/Pthy_HT2_Phase_II)               # Thymus 
         Cadp_HT2_Phase_II <- Aadp_HT2_Phase_II/Vadp
         dAadp_HT2_Phase_II <- Qadp*(Ca_HT2_Phase_II - Cadp_HT2_Phase_II/Padp_HT2_Phase_II)               # Adipose tissue + yellow marrow
         Cmus_HT2_Phase_II <- Amus_HT2_Phase_II/Vmus
         dAmus_HT2_Phase_II <- Qmus*(Ca_HT2_Phase_II - Cmus_HT2_Phase_II/Pmus_HT2_Phase_II)               # Skeletal muscle
         Cbon_HT2_Phase_II <- Abon_HT2_Phase_II/Vbon
         dAbon_HT2_Phase_II <- Qbon*(Ca_HT2_Phase_II - Cbon_HT2_Phase_II/Pbon_HT2_Phase_II)               # Bone + red marrow
         Cgon_HT2_Phase_II <- Agon_HT2_Phase_II/Vgon
         dAgon_HT2_Phase_II <- Qgon*(Ca_HT2_Phase_II - Cgon_HT2_Phase_II/Pgon_HT2_Phase_II)               # Gonads
         Cpan_HT2_Phase_II <- Apan_HT2_Phase_II/Vpan
         dApan_HT2_Phase_II <- Qpan*(Ca_HT2_Phase_II - Cpan_HT2_Phase_II/Ppan_HT2_Phase_II)               # Pancreas
         Cspl_HT2_Phase_II <- Aspl_HT2_Phase_II/Vspl
         dAspl_HT2_Phase_II <- Qspl*(Ca_HT2_Phase_II - Cspl_HT2_Phase_II/Pspl_HT2_Phase_II)               # Spleen
         Crem_HT2_Phase_II <- Arem_HT2_Phase_II/Vrem
         dArem_HT2_Phase_II <- Qrem*(Ca_HT2_Phase_II - Crem_HT2_Phase_II/Prem_HT2_Phase_II)               # Remaining tissues
         
         # Clearing compartments
         
         # Liver
         
         Cport_HT2_Phase_II <- (Qsto*(Csto_HT2_Phase_II/Psto_HT2_Phase_II) +  # Portal blood 
                                  Qint*(Cint_HT2_Phase_II/Pint_HT2_Phase_II) +                            # concentration
                                  Qpan*(Cpan_HT2_Phase_II/Ppan_HT2_Phase_II) + 
                                  Qspl*(Cspl_HT2_Phase_II/Pspl_HT2_Phase_II))/Qport  
         Cliv_HT2_Phase_II <- Aliv_HT2_Phase_II/Vliv  
         dAliv_HT2_Phase_II <- Qliv*Ca_HT2_Phase_II + Qport*Cport_HT2_Phase_II -                       # Liver 
           (Qliv + Qport)*(Cliv_HT2_Phase_II/Pliv_HT2_Phase_II) +      
           CLliv_HT2_HT2_Phase_II*(Cliv_HT2/Pliv_HT2) +
           CLliv_T2_Phase_II_HT2_Phase_II*(Cliv_T2_Phase_II/Pliv_T2_Phase_II) # metabolic conversion of T2 into HT2_Phase_II
         
         dAlivCL_HT2_Phase_II <- 0
         # or: Vmaxliv*(Cliv/Pliv)/(kMliv + Cliv/Pliv)     
         
         # Kidneys
         
         Ckid_HT2_Phase_II <- Akid_HT2_Phase_II/Vkid
         dAkid_HT2_Phase_II <- Qkid*(Ca_HT2_Phase_II - Ckid_HT2_Phase_II/Pkid_HT2_Phase_II) - CLkid_HT2_Phase_II*Ca_HT2_Phase_II        # Kidneys
         dAkidCL_HT2_Phase_II <- CLkid_HT2_Phase_II*Ca_HT2_Phase_II
         
         # Arterial and venous blood
         
         # Note: 5% of Qc flows through non-included tissues (lymph, a.o.) 
         # TODO: check balance Q
         
         dAa_HT2_Phase_II <- ifelse(include_lung==TRUE, 
                                    Qc*(Clun_HT2_Phase_II/Plunb_HT2_Phase_II),
                                    Qc*(Calv_HT2_Phase_II/Pab_HT2_Phase_II)) -                               # Arterial blood
           (QskE + QskU + Qhrt + Qbrn + Qthy + Qadp +
              Qmus + Qbon + Qgon + Qrem + Qsto + Qint + 
              Qpan + Qspl + Qliv + Qkid + 0.05*Qc)*Ca_HT2_Phase_II
         
         dAv_HT2_Phase_II <- Aiv_HT2_Phase_II +                                           # Venous blood
           QskE*CskE_HT2_Phase_II/Pskn_HT2_Phase_II + 
           QskU*CskU_HT2_Phase_II/Pskn_HT2_Phase_II + 
           Qhrt*Chrt_HT2_Phase_II/Phrt_HT2_Phase_II + 
           Qbrn*Cbrn_HT2_Phase_II/Pbrn_HT2_Phase_II + 
           Qthy*Cthy_HT2_Phase_II/Pthy_HT2_Phase_II + 
           Qadp*Cadp_HT2_Phase_II/Padp_HT2_Phase_II + 
           Qmus*Cmus_HT2_Phase_II/Pmus_HT2_Phase_II + 
           Qbon*Cbon_HT2_Phase_II/Pbon_HT2_Phase_II + 
           Qgon*Cgon_HT2_Phase_II/Pgon_HT2_Phase_II + 
           Qrem*Crem_HT2_Phase_II/Prem_HT2_Phase_II + 
           (Qliv + Qport)*(Cliv_HT2_Phase_II/Pliv_HT2_Phase_II) +
           Qkid*Ckid_HT2_Phase_II/Pkid_HT2_Phase_II +
           0.05*Qc*Ca_HT2_Phase_II - Qc*Cv_HT2_Phase_II                    
         
         # Mass balance
         dCumulativeIntake_HT2_Phase_II <- intake_HT2_Phase_II #+ km_sto_T2*AstoL_T2 + km_int_T2*AstoL_T2 + CLliv_T2*(Cliv_T2/Pliv_T2) 
         MassBal_HT2_Phase_II <- (Aalv_HT2_Phase_II + Alun_HT2_Phase_II + Ass_HT2_Phase_II + AskE_HT2_Phase_II + AskU_HT2_Phase_II + AstoL_HT2_Phase_II + 
                                    AintL_HT2_Phase_II + Aint_metab_HT2_Phase_II + Asto_HT2_Phase_II + Aint_HT2_Phase_II + Afec_HT2_Phase_II + Ahrt_HT2_Phase_II + Abrn_HT2_Phase_II + Athy_HT2_Phase_II + Aadp_HT2_Phase_II + Amus_HT2_Phase_II + Abon_HT2_Phase_II + 
                                    Agon_HT2_Phase_II + Arem_HT2_Phase_II + Apan_HT2_Phase_II + Aspl_HT2_Phase_II + Aliv_HT2_Phase_II + AlivCL_HT2_Phase_II + Akid_HT2_Phase_II + AkidCL_HT2_Phase_II + Aa_HT2_Phase_II + Av_HT2_Phase_II)#/CumulativeIntake_HT2_Phase_II 
         
         # Total T-2 and HT-2 intake
         Atotal <- CumulativeIntake_T2 + CumulativeIntake_HT2
         
         # Amount of T-2 and HT-2 in urine includes unchanged and glucuronidated T-2 and HT-2
         Aurine_T2 <- AkidCL_T2 + AkidCL_T2_Phase_II
         Aurine_HT2 <- AkidCL_HT2 + AkidCL_HT2_Phase_II
         # Output
         
         list(c(dAalv_T2, dAlun_T2, dAss_T2, dAskE_T2, dAskU_T2, dAoral_total_T2, dAstoL_T2, dAintL_T2, dAint_metab_T2, dAsto_T2, dAint_T2, dAfec_T2, 
                dAhrt_T2, dAbrn_T2, dAthy_T2, dAadp_T2, dAmus_T2, dAbon_T2, dAgon_T2, dArem_T2, dApan_T2, dAspl_T2, 
                dAliv_T2, dAlivCL_T2, dAkid_T2, dAkidCL_T2, dAa_T2, dAv_T2, dCumulativeIntake_T2,
                
                dAalv_HT2, dAlun_HT2, dAss_HT2, dAskE_HT2, dAskU_HT2, dAoral_total_HT2, dAstoL_HT2, dAintL_HT2, dAint_metab_HT2, dAsto_HT2, dAint_HT2, dAfec_HT2, 
                dAhrt_HT2, dAbrn_HT2, dAthy_HT2, dAadp_HT2, dAmus_HT2, dAbon_HT2, dAgon_HT2, dArem_HT2, dApan_HT2, dAspl_HT2, 
                dAliv_HT2, dAlivCL_HT2, dAkid_HT2, dAkidCL_HT2, dAa_HT2, dAv_HT2, dCumulativeIntake_HT2,
                
                dAalv_T2_Phase_I, dAlun_T2_Phase_I, dAss_T2_Phase_I, dAskE_T2_Phase_I, dAskU_T2_Phase_I, dAoral_total_T2_Phase_I, dAstoL_T2_Phase_I, dAintL_T2_Phase_I, dAint_metab_T2_Phase_I, dAsto_T2_Phase_I, dAint_T2_Phase_I, dAfec_T2_Phase_I, 
                dAhrt_T2_Phase_I, dAbrn_T2_Phase_I, dAthy_T2_Phase_I, dAadp_T2_Phase_I, dAmus_T2_Phase_I, dAbon_T2_Phase_I, dAgon_T2_Phase_I, dArem_T2_Phase_I, dApan_T2_Phase_I, dAspl_T2_Phase_I, 
                dAliv_T2_Phase_I, dAlivCL_T2_Phase_I, dAkid_T2_Phase_I, dAkidCL_T2_Phase_I, dAa_T2_Phase_I, dAv_T2_Phase_I, dCumulativeIntake_T2_Phase_I,
                
                dAalv_HT2_Phase_I, dAlun_HT2_Phase_I, dAss_HT2_Phase_I, dAskE_HT2_Phase_I, dAskU_HT2_Phase_I, dAoral_total_HT2_Phase_I, dAstoL_HT2_Phase_I, dAintL_HT2_Phase_I, dAint_metab_HT2_Phase_I, dAsto_HT2_Phase_I, dAint_HT2_Phase_I, dAfec_HT2_Phase_I, 
                dAhrt_HT2_Phase_I, dAbrn_HT2_Phase_I, dAthy_HT2_Phase_I, dAadp_HT2_Phase_I, dAmus_HT2_Phase_I, dAbon_HT2_Phase_I, dAgon_HT2_Phase_I, dArem_HT2_Phase_I, dApan_HT2_Phase_I, dAspl_HT2_Phase_I, 
                dAliv_HT2_Phase_I, dAlivCL_HT2_Phase_I, dAkid_HT2_Phase_I, dAkidCL_HT2_Phase_I, dAa_HT2_Phase_I, dAv_HT2_Phase_I, dCumulativeIntake_HT2_Phase_I,
                
                dAalv_T2_Phase_II, dAlun_T2_Phase_II, dAss_T2_Phase_II, dAskE_T2_Phase_II, dAskU_T2_Phase_II, dAoral_total_T2_Phase_II, dAstoL_T2_Phase_II, dAintL_T2_Phase_II, dAint_metab_T2_Phase_II, dAsto_T2_Phase_II, dAint_T2_Phase_II, dAfec_T2_Phase_II, 
                dAhrt_T2_Phase_II, dAbrn_T2_Phase_II, dAthy_T2_Phase_II, dAadp_T2_Phase_II, dAmus_T2_Phase_II, dAbon_T2_Phase_II, dAgon_T2_Phase_II, dArem_T2_Phase_II, dApan_T2_Phase_II, dAspl_T2_Phase_II, 
                dAliv_T2_Phase_II, dAlivCL_T2_Phase_II, dAkid_T2_Phase_II, dAkidCL_T2_Phase_II, dAa_T2_Phase_II, dAv_T2_Phase_II, dCumulativeIntake_T2_Phase_II,
                
                dAalv_HT2_Phase_II, dAlun_HT2_Phase_II, dAss_HT2_Phase_II, dAskE_HT2_Phase_II, dAskU_HT2_Phase_II, dAoral_total_HT2_Phase_II, dAstoL_HT2_Phase_II, dAintL_HT2_Phase_II, dAint_metab_HT2_Phase_II, dAsto_HT2_Phase_II, dAint_HT2_Phase_II, dAfec_HT2_Phase_II, 
                dAhrt_HT2_Phase_II, dAbrn_HT2_Phase_II, dAthy_HT2_Phase_II, dAadp_HT2_Phase_II, dAmus_HT2_Phase_II, dAbon_HT2_Phase_II, dAgon_HT2_Phase_II, dArem_HT2_Phase_II, dApan_HT2_Phase_II, dAspl_HT2_Phase_II, 
                dAliv_HT2_Phase_II, dAlivCL_HT2_Phase_II, dAkid_HT2_Phase_II, dAkidCL_HT2_Phase_II, dAa_HT2_Phase_II, dAv_HT2_Phase_II, dCumulativeIntake_HT2_Phase_II),
              c(MassBal_T2=MassBal_T2,Fabs_T2=Fabs_T2,
                Calv_T2=Calv_T2, Css_T2=Css_T2, CskE_T2=CskE_T2, CskU_T2=CskU_T2, 
                Csto_T2=Csto_T2, Cint_T2=Cint_T2, Chrt_T2=Chrt_T2,Cbrn_T2=Cbrn_T2,Cthy_T2=Cthy_T2,
                Cadp_T2=Cadp_T2,Cmus_T2=Cmus_T2,Cbon_T2=Cbon_T2,Cgon_T2=Cgon_T2,Cpan_T2=Cpan_T2,Cspl_T2=Cspl_T2,
                Crem_T2=Crem_T2,Ckid_T2=Ckid_T2,Cliv_T2=Cliv_T2,Ca_T2=Ca_T2,Cv_T2=Cv_T2,
                Aair_T2=Aair_T2,As_T2=As_T2,Aoral_T2=Aoral_T2,Aiv_T2=Aiv_T2, 
                
                MassBal_HT2=MassBal_HT2,Fabs_HT2=Fabs_HT2,
                Calv_HT2=Calv_HT2, Css_HT2=Css_HT2, CskE_HT2=CskE_HT2, CskU_HT2=CskU_HT2, 
                Csto_HT2=Csto_HT2, Cint_HT2=Cint_HT2, Chrt_HT2=Chrt_HT2,Cbrn_HT2=Cbrn_HT2,Cthy_HT2=Cthy_HT2,
                Cadp_HT2=Cadp_HT2,Cmus_HT2=Cmus_HT2,Cbon_HT2=Cbon_HT2,Cgon_HT2=Cgon_HT2,Cpan_HT2=Cpan_HT2,Cspl_HT2=Cspl_HT2,
                Crem_HT2=Crem_HT2,Ckid_HT2=Ckid_HT2,Cliv_HT2=Cliv_HT2,Ca_HT2=Ca_HT2,Cv_HT2=Cv_HT2,
                Aair_HT2=Aair_HT2,As_HT2=As_HT2,Aoral_HT2=Aoral_HT2,Aiv_HT2=Aiv_HT2, 
                
                MassBal_T2_Phase_I=MassBal_T2_Phase_I,Fabs_T2_Phase_I=Fabs_T2_Phase_I,
                Calv_T2_Phase_I=Calv_T2_Phase_I, Css_T2_Phase_I=Css_T2_Phase_I, CskE_T2_Phase_I=CskE_T2_Phase_I, CskU_T2_Phase_I=CskU_T2_Phase_I, 
                Csto_T2_Phase_I=Csto_T2_Phase_I, Cint_T2_Phase_I=Cint_T2_Phase_I, Chrt_T2_Phase_I=Chrt_T2_Phase_I,Cbrn_T2_Phase_I=Cbrn_T2_Phase_I,Cthy_T2_Phase_I=Cthy_T2_Phase_I,
                Cadp_T2_Phase_I=Cadp_T2_Phase_I,Cmus_T2_Phase_I=Cmus_T2_Phase_I,Cbon_T2_Phase_I=Cbon_T2_Phase_I,Cgon_T2_Phase_I=Cgon_T2_Phase_I,Cpan_T2_Phase_I=Cpan_T2_Phase_I,Cspl_T2_Phase_I=Cspl_T2_Phase_I,
                Crem_T2_Phase_I=Crem_T2_Phase_I,Ckid_T2_Phase_I=Ckid_T2_Phase_I,Cliv_T2_Phase_I=Cliv_T2_Phase_I,Ca_T2_Phase_I=Ca_T2_Phase_I,Cv_T2_Phase_I=Cv_T2_Phase_I,
                Aair_T2_Phase_I=Aair_T2_Phase_I,As_T2_Phase_I=As_T2_Phase_I,Aoral_T2_Phase_I=Aoral_T2_Phase_I,Aiv_T2_Phase_I=Aiv_T2_Phase_I,
                
                MassBal_HT2_Phase_I=MassBal_HT2_Phase_I,Fabs_HT2_Phase_I=Fabs_HT2_Phase_I,
                Calv_HT2_Phase_I=Calv_HT2_Phase_I, Css_HT2_Phase_I=Css_HT2_Phase_I, CskE_HT2_Phase_I=CskE_HT2_Phase_I, CskU_HT2_Phase_I=CskU_HT2_Phase_I, 
                Csto_HT2_Phase_I=Csto_HT2_Phase_I, Cint_HT2_Phase_I=Cint_HT2_Phase_I, Chrt_HT2_Phase_I=Chrt_HT2_Phase_I,Cbrn_HT2_Phase_I=Cbrn_HT2_Phase_I,Cthy_HT2_Phase_I=Cthy_HT2_Phase_I,
                Cadp_HT2_Phase_I=Cadp_HT2_Phase_I,Cmus_HT2_Phase_I=Cmus_HT2_Phase_I,Cbon_HT2_Phase_I=Cbon_HT2_Phase_I,Cgon_HT2_Phase_I=Cgon_HT2_Phase_I,Cpan_HT2_Phase_I=Cpan_HT2_Phase_I,Cspl_HT2_Phase_I=Cspl_HT2_Phase_I,
                Crem_HT2_Phase_I=Crem_HT2_Phase_I,Ckid_HT2_Phase_I=Ckid_HT2_Phase_I,Cliv_HT2_Phase_I=Cliv_HT2_Phase_I,Ca_HT2_Phase_I=Ca_HT2_Phase_I,Cv_HT2_Phase_I=Cv_HT2_Phase_I,
                Aair_HT2_Phase_I=Aair_HT2_Phase_I,As_HT2_Phase_I=As_HT2_Phase_I,Aoral_HT2_Phase_I=Aoral_HT2_Phase_I,Aiv_HT2_Phase_I=Aiv_HT2_Phase_I,
                
                MassBal_T2_Phase_II=MassBal_T2_Phase_II,Fabs_T2_Phase_II=Fabs_T2_Phase_II,
                Calv_T2_Phase_II=Calv_T2_Phase_II, Css_T2_Phase_II=Css_T2_Phase_II, CskE_T2_Phase_II=CskE_T2_Phase_II, CskU_T2_Phase_II=CskU_T2_Phase_II, 
                Csto_T2_Phase_II=Csto_T2_Phase_II, Cint_T2_Phase_II=Cint_T2_Phase_II, Chrt_T2_Phase_II=Chrt_T2_Phase_II,Cbrn_T2_Phase_II=Cbrn_T2_Phase_II,Cthy_T2_Phase_II=Cthy_T2_Phase_II,
                Cadp_T2_Phase_II=Cadp_T2_Phase_II,Cmus_T2_Phase_II=Cmus_T2_Phase_II,Cbon_T2_Phase_II=Cbon_T2_Phase_II,Cgon_T2_Phase_II=Cgon_T2_Phase_II,Cpan_T2_Phase_II=Cpan_T2_Phase_II,Cspl_T2_Phase_II=Cspl_T2_Phase_II,
                Crem_T2_Phase_II=Crem_T2_Phase_II,Ckid_T2_Phase_II=Ckid_T2_Phase_II,Cliv_T2_Phase_II=Cliv_T2_Phase_II,Ca_T2_Phase_II=Ca_T2_Phase_II,Cv_T2_Phase_II=Cv_T2_Phase_II,
                Aair_T2_Phase_II=Aair_T2_Phase_II,As_T2_Phase_II=As_T2_Phase_II,Aoral_T2_Phase_II=Aoral_T2_Phase_II,Aiv_T2_Phase_II=Aiv_T2_Phase_II,
                
                MassBal_HT2_Phase_II=MassBal_HT2_Phase_II,Fabs_HT2_Phase_II=Fabs_HT2_Phase_II,
                Calv_HT2_Phase_II=Calv_HT2_Phase_II, Css_HT2_Phase_II=Css_HT2_Phase_II, CskE_HT2_Phase_II=CskE_HT2_Phase_II, CskU_HT2_Phase_II=CskU_HT2_Phase_II, 
                Csto_HT2_Phase_II=Csto_HT2_Phase_II, Cint_HT2_Phase_II=Cint_HT2_Phase_II, Chrt_HT2_Phase_II=Chrt_HT2_Phase_II,Cbrn_HT2_Phase_II=Cbrn_HT2_Phase_II,Cthy_HT2_Phase_II=Cthy_HT2_Phase_II,
                Cadp_HT2_Phase_II=Cadp_HT2_Phase_II,Cmus_HT2_Phase_II=Cmus_HT2_Phase_II,Cbon_HT2_Phase_II=Cbon_HT2_Phase_II,Cgon_HT2_Phase_II=Cgon_HT2_Phase_II,Cpan_HT2_Phase_II=Cpan_HT2_Phase_II,Cspl_HT2_Phase_II=Cspl_HT2_Phase_II,
                Crem_HT2_Phase_II=Crem_HT2_Phase_II,Ckid_HT2_Phase_II=Ckid_HT2_Phase_II,Cliv_HT2_Phase_II=Cliv_HT2_Phase_II,Ca_HT2_Phase_II=Ca_HT2_Phase_II,Cv_HT2_Phase_II=Cv_HT2_Phase_II,
                Aair_HT2_Phase_II=Aair_HT2_Phase_II,As_HT2_Phase_II=As_HT2_Phase_II,Aoral_HT2_Phase_II=Aoral_HT2_Phase_II,Aiv_HT2_Phase_II=Aiv_HT2_Phase_II,
                
                Qp=Qp,Atotal=Atotal,Aurine_T2=Aurine_T2,Aurine_HT2=Aurine_HT2
              ))
       })
}
