f.parms.pbk <- function(phyB = phB, Ptib_T2 = Ptb_T2, Ptib_HT2 = Ptb_HT2, 
                        Ptib_T2_Phase_II = Ptb_T2_Phase_II, Ptib_HT2_Phase_II = Ptb_HT2_Phase_II, 
                        Ptib_T2_Phase_I = Ptb_T2_Phase_I, Ptib_HT2_Phase_I = Ptb_HT2_Phase_I)
{
  # Generates a list of parameter values as input to the PBK model f.pbk. 
  # Physiological parameters are derived from f.physB, partitioning coefficients 
  # from f.Ptisb. 
  
  parms <- c(BSA    = NA,      # 1 - body surface area
             Qc     = NA,      # 2 - cardiac output
             Vlun   = NA,      # 3 - lung tissue volume
             Vhrt   = NA,      # 4 - heart volume
             Vskn   = NA,      # 5 - skin volume
             Vadp   = NA,      # 6 - adipose tissue volume
             Vmus   = NA,      # 7 - muscle volume
             Vbon   = NA,      # 8 - bone volume
             Vbrn   = NA,      # 9 - brain volume
             Vthy   = NA,      # 10 - thymus volume
             Vgon   = NA,      # 11 - gonads volume
             Vkid   = NA,      # 12 - kidney volume
             Vsto   = NA,      # 13 - stomach volume
             Vint   = NA,      # 14 - intestine volume
             Vspl   = NA,      # 15 - spleen volume
             Vpan   = NA,      # 16 - pancreas volume
             Vliv   = NA,      # 17 - liver volume
             Vrem   = NA,      # 18 - remaining tissues volume
             Vb     = NA,      # 19 - blood volume
             Qhrt   = NA,      # 20 - blood flow to heart
             Qskn   = NA,      # 21 - blood flow to skin
             Qadp   = NA,      # 22 - blood flow to adipose tissue
             Qmus   = NA,      # 23 - blood flow to muscle
             Qbon   = NA,      # 24 - blood flow to bone
             Qbrn   = NA,      # 25 - blood flow to brain
             Qthy   = NA,      # 26 - blood flow to thymus
             Qgon   = NA,      # 27 - blood flow to gonads
             Qkid   = NA,      # 28 - blood flow to kidney
             Qsto   = NA,      # 29 - blood flow to stomach
             Qint   = NA,      # 30 - blood flow to intestine
             Qspl   = NA,      # 31 - blood flow to spleen
             Qpan   = NA,      # 32 - blood flow to pancreas
             Qliv   = NA,      # 33 - blood flow to liver
             Qrem   = NA,      # 34 - blood flow to remaining tissues
             Qp     = NA,      # 35 - ventilation rate
             Valv   = NA,      # 36 - alveolar volume
             Plunb_T2  = NA,      # 37 - lung:blood partition coefficient
             Phrt_T2   = NA,      # 38 - heart:blood partition coefficient
             Pskn_T2   = NA,      # 39 - skin:blood partition coefficient
             Padp_T2   = NA,      # 40 - adipose tissue:blood partition coefficient 
             Pmus_T2   = NA,      # 41 - muscle:blood partition coefficient 
             Pbon_T2   = NA,      # 42 - bone:blood partition coefficient 
             Pbrn_T2   = NA,      # 43 - brain:blood partition coefficient
             Pthy_T2   = NA,      # 44 - thymus:blood partition coefficient  
             Pgon_T2   = NA,      # 45 - gonads:blood partition coefficient 
             Pkid_T2   = NA,      # 46 - kidney:blood partition coefficient
             Psto_T2   = NA,      # 47 - stomach:blood partition coefficient
             Pint_T2   = NA,      # 48 - intestine:blood partition coefficient
             Pspl_T2   = NA,      # 49 - spleen:blood partition coefficient 
             Ppan_T2   = NA,      # 50 - pancreas:blood partition coefficient 
             Pliv_T2   = NA,      # 51 - liver:blood partition coefficient 
             Prem_T2   = NA,      # 52 - remaining tissues:blood partition coefficient 
             
             Plunb_HT2  = NA,      # 53 - lung:blood partition coefficient
             Phrt_HT2   = NA,      # 54 - heart:blood partition coefficient
             Pskn_HT2   = NA,      # 55 - skin:blood partition coefficient
             Padp_HT2   = NA,      # 56 - adipose tissue:blood partition coefficient 
             Pmus_HT2   = NA,      # 57 - muscle:blood partition coefficient 
             Pbon_HT2   = NA,      # 58 - bone:blood partition coefficient 
             Pbrn_HT2   = NA,      # 59 - brain:blood partition coefficient
             Pthy_HT2   = NA,      # 60 - thymus:blood partition coefficient  
             Pgon_HT2   = NA,      # 61 - gonads:blood partition coefficient 
             Pkid_HT2   = NA,      # 62 - kidney:blood partition coefficient
             Psto_HT2   = NA,      # 63 - stomach:blood partition coefficient
             Pint_HT2   = NA,      # 64 - intestine:blood partition coefficient
             Pspl_HT2   = NA,      # 65 - spleen:blood partition coefficient 
             Ppan_HT2   = NA,      # 66 - pancreas:blood partition coefficient 
             Pliv_HT2   = NA,      # 67 - liver:blood partition coefficient 
             Prem_HT2   = NA,      # 68 - remaining tissues:blood partition coefficient 
             
             Plunb_T2_Phase_I  = NA,      # 69 - lung:blood partition coefficient
             Phrt_T2_Phase_I   = NA,      # 70 - heart:blood partition coefficient
             Pskn_T2_Phase_I   = NA,      # 71 - skin:blood partition coefficient
             Padp_T2_Phase_I   = NA,      # 72 - adipose tissue:blood partition coefficient 
             Pmus_T2_Phase_I   = NA,      # 73 - muscle:blood partition coefficient 
             Pbon_T2_Phase_I   = NA,      # 74 - bone:blood partition coefficient 
             Pbrn_T2_Phase_I   = NA,      # 75 - brain:blood partition coefficient
             Pthy_T2_Phase_I   = NA,      # 76 - thymus:blood partition coefficient  
             Pgon_T2_Phase_I   = NA,      # 77 - gonads:blood partition coefficient 
             Pkid_T2_Phase_I   = NA,      # 78 - kidney:blood partition coefficient
             Psto_T2_Phase_I   = NA,      # 79 - stomach:blood partition coefficient
             Pint_T2_Phase_I   = NA,      # 80 - intestine:blood partition coefficient
             Pspl_T2_Phase_I   = NA,      # 81 - spleen:blood partition coefficient 
             Ppan_T2_Phase_I   = NA,      # 82 - pancreas:blood partition coefficient 
             Pliv_T2_Phase_I   = NA,      # 83 - liver:blood partition coefficient 
             Prem_T2_Phase_I   = NA,      # 84 - remaining tissues:blood partition coefficient 
             
             Plunb_HT2_Phase_I  = NA,      # 85 - lung:blood partition coefficient
             Phrt_HT2_Phase_I   = NA,      # 86 - heart:blood partition coefficient
             Pskn_HT2_Phase_I   = NA,      # 87 - skin:blood partition coefficient
             Padp_HT2_Phase_I   = NA,      # 88 - adipose tissue:blood partition coefficient 
             Pmus_HT2_Phase_I   = NA,      # 89 - muscle:blood partition coefficient 
             Pbon_HT2_Phase_I   = NA,      # 90 - bone:blood partition coefficient 
             Pbrn_HT2_Phase_I   = NA,      # 91 - brain:blood partition coefficient
             Pthy_HT2_Phase_I   = NA,      # 92 - thymus:blood partition coefficient  
             Pgon_HT2_Phase_I   = NA,      # 93 - gonads:blood partition coefficient 
             Pkid_HT2_Phase_I   = NA,      # 94 - kidney:blood partition coefficient
             Psto_HT2_Phase_I   = NA,      # 95 - stomach:blood partition coefficient
             Pint_HT2_Phase_I   = NA,      # 96 - intestine:blood partition coefficient
             Pspl_HT2_Phase_I   = NA,      # 97 - spleen:blood partition coefficient 
             Ppan_HT2_Phase_I   = NA,      # 98 - pancreas:blood partition coefficient 
             Pliv_HT2_Phase_I   = NA,      # 99 - liver:blood partition coefficient 
             Prem_HT2_Phase_I   = NA,      # 100 - remaining tissues:blood partition coefficient 
             
             Plunb_T2_Phase_II  = NA,      # 101 - lung:blood partition coefficient
             Phrt_T2_Phase_II   = NA,      # 102 - heart:blood partition coefficient
             Pskn_T2_Phase_II   = NA,      # 103 - skin:blood partition coefficient
             Padp_T2_Phase_II   = NA,      # 104 - adipose tissue:blood partition coefficient 
             Pmus_T2_Phase_II   = NA,      # 105 - muscle:blood partition coefficient 
             Pbon_T2_Phase_II   = NA,      # 106 - bone:blood partition coefficient 
             Pbrn_T2_Phase_II   = NA,      # 107 - brain:blood partition coefficient
             Pthy_T2_Phase_II   = NA,      # 108 - thymus:blood partition coefficient  
             Pgon_T2_Phase_II   = NA,      # 109 - gonads:blood partition coefficient 
             Pkid_T2_Phase_II   = NA,      # 110 - kidney:blood partition coefficient
             Psto_T2_Phase_II   = NA,      # 111 - stomach:blood partition coefficient
             Pint_T2_Phase_II   = NA,      # 112 - intestine:blood partition coefficient
             Pspl_T2_Phase_II   = NA,      # 113 - spleen:blood partition coefficient 
             Ppan_T2_Phase_II   = NA,      # 114 - pancreas:blood partition coefficient 
             Pliv_T2_Phase_II   = NA,      # 115 - liver:blood partition coefficient 
             Prem_T2_Phase_II   = NA,      # 116 - remaining tissues:blood partition coefficient 
             
             Plunb_HT2_Phase_II  = NA,      # 117 - lung:blood partition coefficient
             Phrt_HT2_Phase_II   = NA,      # 118 - heart:blood partition coefficient
             Pskn_HT2_Phase_II   = NA,      # 119 - skin:blood partition coefficient
             Padp_HT2_Phase_II   = NA,      # 120 - adipose tissue:blood partition coefficient 
             Pmus_HT2_Phase_II   = NA,      # 121 - muscle:blood partition coefficient 
             Pbon_HT2_Phase_II   = NA,      # 122 - bone:blood partition coefficient 
             Pbrn_HT2_Phase_II   = NA,      # 123 - brain:blood partition coefficient
             Pthy_HT2_Phase_II   = NA,      # 124 - thymus:blood partition coefficient  
             Pgon_HT2_Phase_II   = NA,      # 125 - gonads:blood partition coefficient 
             Pkid_HT2_Phase_II   = NA,      # 126 - kidney:blood partition coefficient
             Psto_HT2_Phase_II   = NA,      # 127 - stomach:blood partition coefficient
             Pint_HT2_Phase_II   = NA,      # 128 - intestine:blood partition coefficient
             Pspl_HT2_Phase_II   = NA,      # 129 - spleen:blood partition coefficient 
             Ppan_HT2_Phase_II   = NA,      # 130 - pancreas:blood partition coefficient 
             Pliv_HT2_Phase_II   = NA,      # 131 - liver:blood partition coefficient 
             Prem_HT2_Phase_II   = NA,      # 132 - remaining tissues:blood partition coefficient 
             
             fss    = 0.,      # 133 - fraction of skin exposed    
             Vss    = 0.1,     # 134 - Matrix volume
             kse    = NA,      # 135 - stomach emptying
             kfec    = NA,     # 136 - stomach emptying
             TL     = 0.,      # 137 - lag time to stomach emptying
             Tse    = NA,      # 138 - stomach transfer time (Not used??)
             
             ka_sto_T2 = 0,      # 139 - stomach absorption rate constant     
             ka_int_T2 = 0,      # 140 - intestinal absorption rate constant     
             km_sto_T2_HT2 = 0,      # 141 - stomach metabolism rate constant     
             km_sto_T2_T2_Phase_II = 0,      # 142 - stomach metabolism rate constant     
             km_sto_T2_T2_Phase_I = 0,      # 143 - stomach metabolism rate constant     
             km_int_T2_HT2 = 0,      # 144 - intestinal metabolism rate constant     
             km_int_T2_T2_Phase_II = 0,      # 145 - intestinal metabolism rate constant     
             km_int_T2_T2_Phase_I = 0,      # 146 - intestinal metabolism rate constant     
             CLliv_T2_HT2 = NA,      # 147 - hepatic clearance per L liver     
             CLliv_T2_T2_Phase_II = NA,      # 148 - hepatic clearance per L liver     
             CLliv_T2_T2_Phase_I = NA,      # 149 - hepatic clearance per L liver     
             CLkid_T2 = NA,      # 150 - renal clearance per L kidney
             PAalv_T2  = NA,      # 151 - permeability area cross product 
             Pluna_T2  = NA,      # 152 - lung:air partition coefficient
             Ks_T2     = 1E-4,    # 153 - skin permeation coefficient
              
             ka_sto_HT2 = 0,      # 154 - stomach absorption rate constant     
             ka_int_HT2 = 0,      # 155 - intestinal absorption rate constant     
             km_sto_HT2_HT2_Phase_II = 0,      # 156 - stomach metabolism rate constant     
             km_sto_HT2_HT2_Phase_I = 0,      # 157 - stomach metabolism rate constant     
             km_int_HT2_HT2_Phase_II = 0,      # 158 - intestinal metabolism rate constant     
             km_int_HT2_HT2_Phase_I = 0,      # 159 - intestinal metabolism rate constant     
             CLliv_HT2_HT2_Phase_II = NA,      # 160 - hepatic clearance per L liver     
             CLliv_HT2_HT2_Phase_I = NA,      # 161 - hepatic clearance per L liver     
             CLkid_HT2 = NA,      # 162 - renal clearance per L kidney
             PAalv_HT2  = NA,      # 163 - permeability area cross product 
             Pluna_HT2  = NA,      # 164 - lung:air partition coefficient
             Ks_HT2     = 1E-4,    # 165 - skin permeation coefficient
             
             ka_sto_T2_Phase_I = 0,      # 166 - stomach absorption rate constant     
             ka_int_T2_Phase_I = 0,      # 167 - intestinal absorption rate constant     
             CLkid_T2_Phase_I = NA,      # 168 - renal clearance per L kidney
             PAalv_T2_Phase_I  = NA,      # 169 - permeability area cross product 
             Pluna_T2_Phase_I  = NA,      # 170 - lung:air partition coefficient
             Ks_T2_Phase_I     = 1E-4,    # 171 - skin permeation coefficient
             
             ka_sto_HT2_Phase_I = 0,      # 172 - stomach absorption rate constant     
             ka_int_HT2_Phase_I = 0,      # 173 - intestinal absorption rate constant     
             CLkid_HT2_Phase_I = NA,      # 174 - renal clearance per L kidney
             PAalv_HT2_Phase_I  = NA,      # 175 - permeability area cross product 
             Pluna_HT2_Phase_I  = NA,      # 176 - lung:air partition coefficient
             Ks_HT2_Phase_I     = 1E-4,    # 177 - skin permeation coefficient
             
             ka_sto_T2_Phase_II = 0,      # 178 - stomach absorption rate constant     
             ka_int_T2_Phase_II = 0,      # 179 - intestinal absorption rate constant     
             km_sto_T2_Phase_II_HT2_Phase_II = 0,      # 180 - stomach metabolism rate constant     
             km_int_T2_Phase_II_HT2_Phase_II = 0,      # 181 - intestinal metabolism rate constant     
             CLliv_T2_Phase_II_HT2_Phase_II = NA,      # 182 - hepatic clearance per L liver     
             CLkid_T2_Phase_II = NA,      # 183 - renal clearance per L kidney
             PAalv_T2_Phase_II  = NA,      # 184 - permeability area cross product 
             Pluna_T2_Phase_II  = NA,      # 185 - lung:air partition coefficient
             Ks_T2_Phase_II     = 1E-4,    # 186 - skin permeation coefficient
             
             ka_sto_HT2_Phase_II = 0,      # 187 - stomach absorption rate constant     
             ka_int_HT2_Phase_II = 0,      # 188 - intestinal absorption rate constant     
             CLkid_HT2_Phase_II = NA,      # 189 - renal clearance per L kidney
             PAalv_HT2_Phase_II  = NA,      # 190 - permeability area cross product 
             Pluna_HT2_Phase_II  = NA,      # 191 - lung:air partition coefficient
             Ks_HT2_Phase_II     = 1E-4    # 192 - skin permeation coefficient
             
)
             
  
  parms[1]     <- phyB[46]      # BSA
  parms[2:34]  <- phyB[4:36]    # Qc, organ volumes, blood flows
  parms[35]    <- phyB[38]      # Qp
  parms[36]    <- phyB[41]      # Valv
  parms[37:52] <- Ptib_T2[1:16]    # Ptis:blood
  parms[53:68] <- Ptib_HT2[1:16]    # Ptis:blood
  parms[69:84] <- Ptib_T2_Phase_I[1:16]    # Ptis:blood
  parms[85:100] <- Ptib_HT2_Phase_I[1:16]   # Ptis:blood
  parms[101:116] <- Ptib_T2_Phase_II[1:16]    # Ptis:blood
  parms[117:132] <- Ptib_HT2_Phase_II[1:16]    # Ptis:blood
  
  return(parms)
}