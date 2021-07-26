#' Vectorized version of Janiczek and DeYoung (1987) for R
#' 
#' @param IY Input year.
#' @param IM Input month
#' @param ID Input day
#' @param LO Longitude
#' @param FINIT LATITUDE
#' @param ZZ Which timezone (0 = GMT, 1 = standard zone time, 2 = local mean time)
#' @param SK Sky condition. (1 = Sun or moon visible, sky less than 70 percentovercast, 2 = sun or moon obscured by thin clouds, 3 = sun or moon obscured by average clouds, 10 = sun/moon obscured by dark stratus clouds (rare))
#' @param HR Time based on a 24 hour clock, as a numeric vector (e.g. 820, 1330) 
#' @param full.output If true, returns sun and moon separately.
#' @return Illuminance (lux) at the Earth's surface.

illumR <- illumR <- function(IY, IM, ID, LO, FINIT, ZZ, SK, HR, full.output = FALSE) {
  
  # Function to convert clock time to decimal time.
  DEG <- function(x) {
    as.integer(x) + ((x-as.integer(x))*100/60)
  }
  
  # Set up parameters-------------------------------------------------------------------------------
  RD <- 57.29577951
  DR <- 1/RD
  AA <- c(-0.01454, -0.10453, -0.20791, 0.00233)
  CE <- 0.91775
  SE <- 0.39715
  
  HH <- HR
  HINIT <- HH
  
  FF <- FINIT
  CC <- 360
  LI <- abs(LO)
  FO <- FF
  FF <- FF * DR
  SI <- sin(FF)
  CI <- cos(FF)
  JJ <- 367*IY-as.integer(7*(IY+as.integer((IM+9)/12))/4) + as.integer(275*IM/9) + ID - 730531
  
  ZT <- ZZ
  DT <- 0
  
  # Deal with timezones-----------------------------------------------------------------------------
  if(ZZ == 0) {
    DT <- -LO/360
  } else if(ZZ == 1) {
    DT <- (LI-15*as.integer((LI+7.5)/15))/CC*sign(-LO)
  }
  
  EE <- DEG(HR/100)/24 - DT - LO/360
  DD <- JJ-0.5+EE
  NN <- 1
  
  # SUN --------------------------------------------------------------------------------------------
  TT <- 280.46 + 0.98565 * DD
  TT <- TT - as.integer(TT/360) *360
  
  if(TT < 0) {
    TT <- TT + 360
  }
  
  GG <- (357.5 + 0.98560 * DD) * DR
  LS <- (TT + 1.91 * sin(GG)) * DR
  T_0 <- LS #
  AS <- atan(CE * tan(LS)) * RD
  YY <- cos(LS)
  
  if(YY < 0) {
    AS <- AS + 180
  }
  
  SD <- SE * sin(LS)
  DS <- asin(SD)
  TT <- TT - 180
  #-------------------------------------------------------------------------------------------------
  
  TT <- TT + 360 * EE + LO
  HH <- TT - AS
  
  # ALTAZ-------------------------------------------------------------------------------------------
  CD <- cos(DS)
  CS <- cos(HH * DR)
  QQ <- SD * CI - CD * SI * CS
  PP <- -CD * sin(HH * DR)
  AZ <- atan(PP/QQ) * RD

  if(QQ < 0) {
    AZ <- AZ + 180
  }
  
  if(AZ < 0) {
    AZ <- AZ + 360
  }
  
  AZ <- AZ + 0.5
  HH <- asin(SD * SI + CD * CI *CS) * RD
  #-------------------------------------------------------------------------------------------------
  ZZ <- HH * DR
  HH <- HH - 0.95 * (NN-1) * cos(HH*DR)
  
  # REFR--------------------------------------------------------------------------------------------
  HA <- HH
  
  if(HH < (-5/6)) {
    # Do nothing
  } else {
    HA <- HH + 1/(tan((HH + 8.6/(HH+4.42))*DR))/60
  }
  # ------------------------------------------------------------------------------------------------
  
  # ATMOS-------------------------------------------------------------------------------------------
  UU <- sin(HA*DR)
  XX <- 753.66156
  SS <- asin(XX * cos(HA*DR)/(XX+1))
  MM <- XX*(cos(SS)-UU)+cos(SS)
  MM <- exp(-0.21*MM)*UU+0.0289*exp(-0.042*MM)*(1.0+(HA+90.0)*UU/57.29577951)
  #-------------------------------------------------------------------------------------------------
  
  HA <- sign(HA)*(abs(HA)+0.5)
  
  IS <- 133775.0 * MM/SK
  IAZ <- AZ
  
  
  # MOON--------------------------------------------------------------------------------------------
  VV <- 218.32 + 13.1764*DD
  VV <- VV - as.integer(VV/360)*360
  if(VV < 0) {
    VV <- VV + 360
  }
  
  YY <- (134.96 + 13.06499 * DD) * DR
  OO <- (93.27 + 13.22935 * DD) * DR
  WW <- (235.7 + 24.38150 * DD) * DR
  SB <- sin(YY)
  CB <- cos(YY)
  XX <- sin(OO)
  SS <- cos(OO)
  SD <- sin(WW)
  CD <- cos(WW)
  VV <- (VV + (6.29-1.27*CD+0.43*CB)*SB+(0.66 + 1.27*CB)*SD-0.19*sin(GG)-0.23*XX*SS)*DR
  YY <- ((5.13-0.17*CD)*XX + (0.56*SB+0.17*SD)*SS)*DR
  SV <- sin(VV)
  SB <- sin(YY)
  CB <- cos(YY)
  QQ <- CB*cos(VV)
  PP <- CE*SV*CB-SE*SB
  SD <- SE*SV*CB+CE*SB
  AS <- atan(PP/QQ)*RD
  if(QQ < 0) {
    AS <- AS+180
  }
  DS <- asin(SD)
  
  # ------------------------------------------------------------------------------------------------
  
  HH <- TT - AS
  
  # ALTAZ-------------------------------------------------------------------------------------------
  CD <- cos(DS)
  CS <- cos(HH * DR)
  QQ <- SD * CI - CD * SI * CS
  PP <- -CD * sin(HH * DR)
  AZ <- atan(PP/QQ) * RD
  
  if(QQ < 0) {
    AZ <- AZ + 180
  }
  
  if(AZ < 0) {
    AZ <- AZ + 360
  }
  
  AZ <- AZ + 0.5
  HH <- asin(SD * SI + CD * CI *CS) * RD
  #-------------------------------------------------------------------------------------------------
  
  ZZ <- HH * DR
  HH <- HH - 0.95 * (NN-1) * cos(HH*DR)
  
  # REFR--------------------------------------------------------------------------------------------
  
  HA <- HH
  
  if(HH < (-5/6)) {
    # Do nothing
  } else {
    HA <- HH + 1/(tan((HH + 8.6/(HH+4.42))*DR))/60
  }
  #-------------------------------------------------------------------------------------------------
  
  # ATMOS-------------------------------------------------------------------------------------------
  UU <- sin(HA*DR)
  XX <- 753.66156
  SS <- asin(XX * cos(HA*DR)/(XX+1))
  MM <- XX*(cos(SS)-UU)+cos(SS)
  MM <- exp(-0.21*MM)*UU+0.0289*exp(-0.042*MM)*(1.0+(HA+90.0)*UU/57.29577951)
  #-------------------------------------------------------------------------------------------------
  
  HA <- sign(HA)*(abs(HA)+0.5)
  
  EE <- acos(cos(VV-LS) * CB)
  PP <- 0.892*exp(-3.343/((tan(EE/2.0))^0.632))+0.0344*(sin(EE)-EE*cos(EE))
  PP <- 0.418*PP/(1.0-0.005*cos(EE)-0.03*sin(ZZ))
  IL <- PP*MM/SK
  ISUN <- IS
  IMOON <- IL
  IS <- IS+IL+0.0005/SK
  IAZ <- AZ
  IHA <- HA
  
  IHA <- 50 * (1-cos(EE)) + 0.5
  
  if(full.output) {
    IS <- data.frame(SUN_ILL = ISUN, MOON_ILL = IMOON, AS = AS, AZ = AZ, CB = CB, CE = CE, DR = DR, EE = EE, GG = GG, HA = HA, HH = HH, HINIT = HINIT, HR = HR, IAZ = IAZ, ID = ID, IHA = IHA, IL = IL, IM = IM, IS = IS,
                     IY = IY, JJ = JJ, LI = LI, LO = LO, LS = LS, MM = MM, NN = NN, OO = OO, PP = PP, QQ = QQ, RD = RD, SB = SB, SD = SD, SE = SE, SI = SI, SK = SK, SS = SS,
                     SV = SV, TT = TT, T_0 = T_0, UU = UU, VV = VV, WW = WW, XX = XX, YY = YY, ZT = ZT, ZZ = ZZ)
  }
  
  return(IS)
}
