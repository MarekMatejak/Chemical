within Chemical;
package Substances "Definitions of substances"
  package IdealGasesMSL

  record Ag "Ag(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ag);
   annotation (preferredView = "info");
  end Ag;

  record Agplus "Agplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Agplus,
    z=1);
   annotation (preferredView = "info");
  end Agplus;

  record Agminus "Agminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Agminus,
    z=-1);
   annotation (preferredView = "info");
  end Agminus;

  record Air "Air(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Air);
   annotation (preferredView = "info");
  end Air;

  record AL "AL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL);
   annotation (preferredView = "info");
  end AL;

  record ALplus "ALplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALplus,
    z=1);
   annotation (preferredView = "info");
  end ALplus;

  record ALminus "ALminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALminus,
    z=-1);
   annotation (preferredView = "info");
  end ALminus;

  record ALBr "ALBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALBr);
   annotation (preferredView = "info");
  end ALBr;

  record ALBr2 "ALBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALBr2);
   annotation (preferredView = "info");
  end ALBr2;

  record ALBr3 "ALBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALBr3);
   annotation (preferredView = "info");
  end ALBr3;

  record ALC "ALC(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALC);
   annotation (preferredView = "info");
  end ALC;

  record ALC2 "ALC2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALC2);
   annotation (preferredView = "info");
  end ALC2;

  record ALCL "ALCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALCL);
   annotation (preferredView = "info");
  end ALCL;

  record ALCLplus "ALCLplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALCLplus,
    z=1);
   annotation (preferredView = "info");
  end ALCLplus;

  record ALCL2 "ALCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALCL2);
   annotation (preferredView = "info");
  end ALCL2;

  record ALCL3 "ALCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALCL3);
   annotation (preferredView = "info");
  end ALCL3;

  record ALF "ALF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF);
   annotation (preferredView = "info");
  end ALF;

  record ALFplus "ALFplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALFplus,
    z=1);
   annotation (preferredView = "info");
  end ALFplus;

  record ALFCL "ALFCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALFCL);
   annotation (preferredView = "info");
  end ALFCL;

  record ALFCL2 "ALFCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALFCL2);
   annotation (preferredView = "info");
  end ALFCL2;

  record ALF2 "ALF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF2);
   annotation (preferredView = "info");
  end ALF2;

  record ALF2minus "ALF2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF2minus,
    z=-1);
   annotation (preferredView = "info");
  end ALF2minus;

  record ALF2CL "ALF2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF2CL);
   annotation (preferredView = "info");
  end ALF2CL;

  record ALF3 "ALF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF3);
   annotation (preferredView = "info");
  end ALF3;

  record ALF4minus "ALF4minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF4minus,
    z=-1);
   annotation (preferredView = "info");
  end ALF4minus;

  record ALH "ALH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH);
   annotation (preferredView = "info");
  end ALH;

  record ALHCL "ALHCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHCL);
   annotation (preferredView = "info");
  end ALHCL;

  record ALHCL2 "ALHCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHCL2);
   annotation (preferredView = "info");
  end ALHCL2;

  record ALHF "ALHF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHF);
   annotation (preferredView = "info");
  end ALHF;

  record ALHFCL "ALHFCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHFCL);
   annotation (preferredView = "info");
  end ALHFCL;

  record ALHF2 "ALHF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHF2);
   annotation (preferredView = "info");
  end ALHF2;

  record ALH2 "ALH2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH2);
   annotation (preferredView = "info");
  end ALH2;

  record ALH2CL "ALH2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH2CL);
   annotation (preferredView = "info");
  end ALH2CL;

  record ALH2F "ALH2F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH2F);
   annotation (preferredView = "info");
  end ALH2F;

  record ALH3 "ALH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH3);
   annotation (preferredView = "info");
  end ALH3;

  record ALI "ALI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALI);
   annotation (preferredView = "info");
  end ALI;

  record ALI2 "ALI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALI2);
   annotation (preferredView = "info");
  end ALI2;

  record ALI3 "ALI3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALI3);
   annotation (preferredView = "info");
  end ALI3;

  record ALN "ALN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALN);
   annotation (preferredView = "info");
  end ALN;

  record ALO "ALO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALO);
   annotation (preferredView = "info");
  end ALO;

  record ALOplus "ALOplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOplus,
    z=1);
   annotation (preferredView = "info");
  end ALOplus;

  record ALOminus "ALOminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOminus,
    z=-1);
   annotation (preferredView = "info");
  end ALOminus;

  record ALOCL "ALOCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOCL);
   annotation (preferredView = "info");
  end ALOCL;

  record ALOCL2 "ALOCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOCL2);
   annotation (preferredView = "info");
  end ALOCL2;

  record ALOF "ALOF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOF);
   annotation (preferredView = "info");
  end ALOF;

  record ALOF2 "ALOF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOF2);
   annotation (preferredView = "info");
  end ALOF2;

  record ALOF2minus "ALOF2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOF2minus,
    z=-1);
   annotation (preferredView = "info");
  end ALOF2minus;

  record ALOH "ALOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOH);
   annotation (preferredView = "info");
  end ALOH;

  record ALOHCL "ALOHCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOHCL);
   annotation (preferredView = "info");
  end ALOHCL;

  record ALOHCL2 "ALOHCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOHCL2);
   annotation (preferredView = "info");
  end ALOHCL2;

  record ALOHF "ALOHF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOHF);
   annotation (preferredView = "info");
  end ALOHF;

  record ALOHF2 "ALOHF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOHF2);
   annotation (preferredView = "info");
  end ALOHF2;

  record ALO2 "ALO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALO2);
   annotation (preferredView = "info");
  end ALO2;

  record ALO2minus "ALO2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALO2minus,
    z=-1);
   annotation (preferredView = "info");
  end ALO2minus;

  record AL_OH_2 "AL_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL_OH_2);
   annotation (preferredView = "info");
  end AL_OH_2;

  record AL_OH_2CL "AL_OH_2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL_OH_2CL);
   annotation (preferredView = "info");
  end AL_OH_2CL;

  record AL_OH_2F "AL_OH_2F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL_OH_2F);
   annotation (preferredView = "info");
  end AL_OH_2F;

  record AL_OH_3 "AL_OH_3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL_OH_3);
   annotation (preferredView = "info");
  end AL_OH_3;

  record ALS "ALS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALS);
   annotation (preferredView = "info");
  end ALS;

  record ALS2 "ALS2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ALS2);
   annotation (preferredView = "info");
  end ALS2;

  record AL2 "AL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2);
   annotation (preferredView = "info");
  end AL2;

  record AL2Br6 "AL2Br6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2Br6);
   annotation (preferredView = "info");
  end AL2Br6;

  record AL2C2 "AL2C2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2C2);
   annotation (preferredView = "info");
  end AL2C2;

  record AL2CL6 "AL2CL6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2CL6);
   annotation (preferredView = "info");
  end AL2CL6;

  record AL2F6 "AL2F6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2F6);
   annotation (preferredView = "info");
  end AL2F6;

  record AL2I6 "AL2I6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2I6);
   annotation (preferredView = "info");
  end AL2I6;

  record AL2O "AL2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2O);
   annotation (preferredView = "info");
  end AL2O;

  record AL2Oplus "AL2Oplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2Oplus,
    z=1);
   annotation (preferredView = "info");
  end AL2Oplus;

  record AL2O2 "AL2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2O2);
   annotation (preferredView = "info");
  end AL2O2;

  record AL2O2plus "AL2O2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2O2plus,
    z=1);
   annotation (preferredView = "info");
  end AL2O2plus;

  record AL2O3 "AL2O3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2O3);
   annotation (preferredView = "info");
  end AL2O3;

  record AL2S "AL2S(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2S);
   annotation (preferredView = "info");
  end AL2S;

  record AL2S2 "AL2S2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2S2);
   annotation (preferredView = "info");
  end AL2S2;

  record Ar "Ar(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ar);
   annotation (preferredView = "info");
  end Ar;

  record Arplus "Arplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Arplus,
    z=1);
   annotation (preferredView = "info");
  end Arplus;

  record B "B(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B);
   annotation (preferredView = "info");
  end B;

  record Bplus "Bplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Bplus,
    z=1);
   annotation (preferredView = "info");
  end Bplus;

  record Bminus "Bminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Bminus,
    z=-1);
   annotation (preferredView = "info");
  end Bminus;

  record BBr "BBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BBr);
   annotation (preferredView = "info");
  end BBr;

  record BBr2 "BBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BBr2);
   annotation (preferredView = "info");
  end BBr2;

  record BBr3 "BBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BBr3);
   annotation (preferredView = "info");
  end BBr3;

  record BC "BC(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BC);
   annotation (preferredView = "info");
  end BC;

  record BC2 "BC2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BC2);
   annotation (preferredView = "info");
  end BC2;

  record BCL "BCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL);
   annotation (preferredView = "info");
  end BCL;

  record BCLplus "BCLplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BCLplus,
    z=1);
   annotation (preferredView = "info");
  end BCLplus;

  record BCLOH "BCLOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BCLOH);
   annotation (preferredView = "info");
  end BCLOH;

  record BCL_OH_2 "BCL_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL_OH_2);
   annotation (preferredView = "info");
  end BCL_OH_2;

  record BCL2 "BCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL2);
   annotation (preferredView = "info");
  end BCL2;

  record BCL2plus "BCL2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL2plus,
    z=1);
   annotation (preferredView = "info");
  end BCL2plus;

  record BCL2OH "BCL2OH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL2OH);
   annotation (preferredView = "info");
  end BCL2OH;

  record BF "BF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BF);
   annotation (preferredView = "info");
  end BF;

  record BFCL "BFCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BFCL);
   annotation (preferredView = "info");
  end BFCL;

  record BFCL2 "BFCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BFCL2);
   annotation (preferredView = "info");
  end BFCL2;

  record BFOH "BFOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BFOH);
   annotation (preferredView = "info");
  end BFOH;

  record BF_OH_2 "BF_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BF_OH_2);
   annotation (preferredView = "info");
  end BF_OH_2;

  record BF2 "BF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2);
   annotation (preferredView = "info");
  end BF2;

  record BF2plus "BF2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2plus,
    z=1);
   annotation (preferredView = "info");
  end BF2plus;

  record BF2minus "BF2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2minus,
    z=-1);
   annotation (preferredView = "info");
  end BF2minus;

  record BF2CL "BF2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2CL);
   annotation (preferredView = "info");
  end BF2CL;

  record BF2OH "BF2OH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2OH);
   annotation (preferredView = "info");
  end BF2OH;

  record BF3 "BF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BF3);
   annotation (preferredView = "info");
  end BF3;

  record BF4minus "BF4minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BF4minus,
    z=-1);
   annotation (preferredView = "info");
  end BF4minus;

  record BH "BH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BH);
   annotation (preferredView = "info");
  end BH;

  record BHCL "BHCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BHCL);
   annotation (preferredView = "info");
  end BHCL;

  record BHCL2 "BHCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BHCL2);
   annotation (preferredView = "info");
  end BHCL2;

  record BHF "BHF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BHF);
   annotation (preferredView = "info");
  end BHF;

  record BHFCL "BHFCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BHFCL);
   annotation (preferredView = "info");
  end BHFCL;

  record BHF2 "BHF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BHF2);
   annotation (preferredView = "info");
  end BHF2;

  record BH2 "BH2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BH2);
   annotation (preferredView = "info");
  end BH2;

  record BH2CL "BH2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BH2CL);
   annotation (preferredView = "info");
  end BH2CL;

  record BH2F "BH2F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BH2F);
   annotation (preferredView = "info");
  end BH2F;

  record BH3 "BH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BH3);
   annotation (preferredView = "info");
  end BH3;

  record BH3NH3 "BH3NH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BH3NH3);
   annotation (preferredView = "info");
  end BH3NH3;

  record BH4 "BH4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BH4);
   annotation (preferredView = "info");
  end BH4;

  record BI "BI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BI);
   annotation (preferredView = "info");
  end BI;

  record BI2 "BI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BI2);
   annotation (preferredView = "info");
  end BI2;

  record BI3 "BI3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BI3);
   annotation (preferredView = "info");
  end BI3;

  record BN "BN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BN);
   annotation (preferredView = "info");
  end BN;

  record BO "BO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BO);
   annotation (preferredView = "info");
  end BO;

  record BOminus "BOminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BOminus,
    z=-1);
   annotation (preferredView = "info");
  end BOminus;

  record BOCL "BOCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BOCL);
   annotation (preferredView = "info");
  end BOCL;

  record BOCL2 "BOCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BOCL2);
   annotation (preferredView = "info");
  end BOCL2;

  record BOF "BOF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BOF);
   annotation (preferredView = "info");
  end BOF;

  record BOF2 "BOF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BOF2);
   annotation (preferredView = "info");
  end BOF2;

  record BOH "BOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BOH);
   annotation (preferredView = "info");
  end BOH;

  record BO2 "BO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BO2);
   annotation (preferredView = "info");
  end BO2;

  record BO2minus "BO2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BO2minus,
    z=-1);
   annotation (preferredView = "info");
  end BO2minus;

  record B_OH_2 "B_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B_OH_2);
   annotation (preferredView = "info");
  end B_OH_2;

  record BS "BS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BS);
   annotation (preferredView = "info");
  end BS;

  record BS2 "BS2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BS2);
   annotation (preferredView = "info");
  end BS2;

  record B2 "B2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2);
   annotation (preferredView = "info");
  end B2;

  record B2C "B2C(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2C);
   annotation (preferredView = "info");
  end B2C;

  record B2CL4 "B2CL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2CL4);
   annotation (preferredView = "info");
  end B2CL4;

  record B2F4 "B2F4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2F4);
   annotation (preferredView = "info");
  end B2F4;

  record B2H "B2H(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H);
   annotation (preferredView = "info");
  end B2H;

  record B2H2 "B2H2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H2);
   annotation (preferredView = "info");
  end B2H2;

  record B2H3 "B2H3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H3);
   annotation (preferredView = "info");
  end B2H3;

  record B2H3_db "B2H3_db(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H3_db);
   annotation (preferredView = "info");
  end B2H3_db;

  record B2H4 "B2H4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H4);
   annotation (preferredView = "info");
  end B2H4;

  record B2H4_db "B2H4_db(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H4_db);
   annotation (preferredView = "info");
  end B2H4_db;

  record B2H5 "B2H5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H5);
   annotation (preferredView = "info");
  end B2H5;

  record B2H5_db "B2H5_db(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H5_db);
   annotation (preferredView = "info");
  end B2H5_db;

  record B2H6 "B2H6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H6);
   annotation (preferredView = "info");
  end B2H6;

  record B2O "B2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2O);
   annotation (preferredView = "info");
  end B2O;

  record B2O2 "B2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2O2);
   annotation (preferredView = "info");
  end B2O2;

  record B2O3 "B2O3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2O3);
   annotation (preferredView = "info");
  end B2O3;

  record B2_OH_4 "B2_OH_4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2_OH_4);
   annotation (preferredView = "info");
  end B2_OH_4;

  record B2S "B2S(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2S);
   annotation (preferredView = "info");
  end B2S;

  record B2S2 "B2S2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2S2);
   annotation (preferredView = "info");
  end B2S2;

  record B2S3 "B2S3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B2S3);
   annotation (preferredView = "info");
  end B2S3;

  record B3H7_C2v "B3H7_C2v(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B3H7_C2v);
   annotation (preferredView = "info");
  end B3H7_C2v;

  record B3H7_Cs "B3H7_Cs(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B3H7_Cs);
   annotation (preferredView = "info");
  end B3H7_Cs;

  record B3H9 "B3H9(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B3H9);
   annotation (preferredView = "info");
  end B3H9;

  record B3N3H6 "B3N3H6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B3N3H6);
   annotation (preferredView = "info");
  end B3N3H6;

  record B3O3CL3 "B3O3CL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B3O3CL3);
   annotation (preferredView = "info");
  end B3O3CL3;

  record B3O3FCL2 "B3O3FCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B3O3FCL2);
   annotation (preferredView = "info");
  end B3O3FCL2;

  record B3O3F2CL "B3O3F2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B3O3F2CL);
   annotation (preferredView = "info");
  end B3O3F2CL;

  record B3O3F3 "B3O3F3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B3O3F3);
   annotation (preferredView = "info");
  end B3O3F3;

  record B4H4 "B4H4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B4H4);
   annotation (preferredView = "info");
  end B4H4;

  record B4H10 "B4H10(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B4H10);
   annotation (preferredView = "info");
  end B4H10;

  record B4H12 "B4H12(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B4H12);
   annotation (preferredView = "info");
  end B4H12;

  record B5H9 "B5H9(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.B5H9);
   annotation (preferredView = "info");
  end B5H9;

  record Ba "Ba(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ba);
   annotation (preferredView = "info");
  end Ba;

  record Baplus "Baplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Baplus,
    z=1);
   annotation (preferredView = "info");
  end Baplus;

  record BaBr "BaBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaBr);
   annotation (preferredView = "info");
  end BaBr;

  record BaBr2 "BaBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaBr2);
   annotation (preferredView = "info");
  end BaBr2;

  record BaCL "BaCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaCL);
   annotation (preferredView = "info");
  end BaCL;

  record BaCLplus "BaCLplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaCLplus,
    z=1);
   annotation (preferredView = "info");
  end BaCLplus;

  record BaCL2 "BaCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaCL2);
   annotation (preferredView = "info");
  end BaCL2;

  record BaF "BaF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaF);
   annotation (preferredView = "info");
  end BaF;

  record BaFplus "BaFplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaFplus,
    z=1);
   annotation (preferredView = "info");
  end BaFplus;

  record BaF2 "BaF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaF2);
   annotation (preferredView = "info");
  end BaF2;

  record BaH "BaH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaH);
   annotation (preferredView = "info");
  end BaH;

  record BaI "BaI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaI);
   annotation (preferredView = "info");
  end BaI;

  record BaI2 "BaI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaI2);
   annotation (preferredView = "info");
  end BaI2;

  record BaO "BaO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaO);
   annotation (preferredView = "info");
  end BaO;

  record BaOplus "BaOplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaOplus,
    z=1);
   annotation (preferredView = "info");
  end BaOplus;

  record BaOH "BaOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaOH);
   annotation (preferredView = "info");
  end BaOH;

  record BaOHplus "BaOHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaOHplus,
    z=1);
   annotation (preferredView = "info");
  end BaOHplus;

  record Ba_OH_2 "Ba_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ba_OH_2);
   annotation (preferredView = "info");
  end Ba_OH_2;

  record BaS "BaS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BaS);
   annotation (preferredView = "info");
  end BaS;

  record Ba2 "Ba2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ba2);
   annotation (preferredView = "info");
  end Ba2;

  record Be "Be(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be);
   annotation (preferredView = "info");
  end Be;

  record Beplus "Beplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Beplus,
    z=1);
   annotation (preferredView = "info");
  end Beplus;

  record Beplusplus "Beplusplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Beplusplus,
    z=1);
   annotation (preferredView = "info");
  end Beplusplus;

  record BeBr "BeBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeBr);
   annotation (preferredView = "info");
  end BeBr;

  record BeBr2 "BeBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeBr2);
   annotation (preferredView = "info");
  end BeBr2;

  record BeCL "BeCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeCL);
   annotation (preferredView = "info");
  end BeCL;

  record BeCL2 "BeCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeCL2);
   annotation (preferredView = "info");
  end BeCL2;

  record BeF "BeF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeF);
   annotation (preferredView = "info");
  end BeF;

  record BeF2 "BeF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeF2);
   annotation (preferredView = "info");
  end BeF2;

  record BeH "BeH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeH);
   annotation (preferredView = "info");
  end BeH;

  record BeHplus "BeHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeHplus,
    z=1);
   annotation (preferredView = "info");
  end BeHplus;

  record BeH2 "BeH2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeH2);
   annotation (preferredView = "info");
  end BeH2;

  record BeI "BeI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeI);
   annotation (preferredView = "info");
  end BeI;

  record BeI2 "BeI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeI2);
   annotation (preferredView = "info");
  end BeI2;

  record BeN "BeN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeN);
   annotation (preferredView = "info");
  end BeN;

  record BeO "BeO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeO);
   annotation (preferredView = "info");
  end BeO;

  record BeOH "BeOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeOH);
   annotation (preferredView = "info");
  end BeOH;

  record BeOHplus "BeOHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeOHplus,
    z=1);
   annotation (preferredView = "info");
  end BeOHplus;

  record Be_OH_2 "Be_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be_OH_2);
   annotation (preferredView = "info");
  end Be_OH_2;

  record BeS "BeS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BeS);
   annotation (preferredView = "info");
  end BeS;

  record Be2 "Be2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2);
   annotation (preferredView = "info");
  end Be2;

  record Be2CL4 "Be2CL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2CL4);
   annotation (preferredView = "info");
  end Be2CL4;

  record Be2F4 "Be2F4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2F4);
   annotation (preferredView = "info");
  end Be2F4;

  record Be2O "Be2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2O);
   annotation (preferredView = "info");
  end Be2O;

  record Be2OF2 "Be2OF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2OF2);
   annotation (preferredView = "info");
  end Be2OF2;

  record Be2O2 "Be2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2O2);
   annotation (preferredView = "info");
  end Be2O2;

  record Be3O3 "Be3O3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be3O3);
   annotation (preferredView = "info");
  end Be3O3;

  record Be4O4 "Be4O4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Be4O4);
   annotation (preferredView = "info");
  end Be4O4;

  record Br "Br(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Br);
   annotation (preferredView = "info");
  end Br;

  record Brplus "Brplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Brplus,
    z=1);
   annotation (preferredView = "info");
  end Brplus;

  record Brminus "Brminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Brminus,
    z=-1);
   annotation (preferredView = "info");
  end Brminus;

  record BrCL "BrCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BrCL);
   annotation (preferredView = "info");
  end BrCL;

  record BrF "BrF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BrF);
   annotation (preferredView = "info");
  end BrF;

  record BrF3 "BrF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BrF3);
   annotation (preferredView = "info");
  end BrF3;

  record BrF5 "BrF5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BrF5);
   annotation (preferredView = "info");
  end BrF5;

  record BrO "BrO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BrO);
   annotation (preferredView = "info");
  end BrO;

  record OBrO "OBrO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.OBrO);
   annotation (preferredView = "info");
  end OBrO;

  record BrOO "BrOO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BrOO);
   annotation (preferredView = "info");
  end BrOO;

  record BrO3 "BrO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BrO3);
   annotation (preferredView = "info");
  end BrO3;

  record Br2 "Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Br2);
   annotation (preferredView = "info");
  end Br2;

  record BrBrO "BrBrO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BrBrO);
   annotation (preferredView = "info");
  end BrBrO;

  record BrOBr "BrOBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.BrOBr);
   annotation (preferredView = "info");
  end BrOBr;

  record C "C(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C);
   annotation (preferredView = "info");
  end C;

  record Cplus "Cplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cplus,
    z=1);
   annotation (preferredView = "info");
  end Cplus;

  record Cminus "Cminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cminus,
    z=-1);
   annotation (preferredView = "info");
  end Cminus;

  record CBr "CBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CBr);
   annotation (preferredView = "info");
  end CBr;

  record CBr2 "CBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CBr2);
   annotation (preferredView = "info");
  end CBr2;

  record CBr3 "CBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CBr3);
   annotation (preferredView = "info");
  end CBr3;

  record CBr4 "CBr4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CBr4);
   annotation (preferredView = "info");
  end CBr4;

  record CCL "CCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL);
   annotation (preferredView = "info");
  end CCL;

  record CCL2 "CCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL2);
   annotation (preferredView = "info");
  end CCL2;

  record CCL2Br2 "CCL2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL2Br2);
   annotation (preferredView = "info");
  end CCL2Br2;

  record CCL3 "CCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL3);
   annotation (preferredView = "info");
  end CCL3;

  record CCL3Br "CCL3Br(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL3Br);
   annotation (preferredView = "info");
  end CCL3Br;

  record CCL4 "CCL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL4);
   annotation (preferredView = "info");
  end CCL4;

  record CF "CF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF);
   annotation (preferredView = "info");
  end CF;

  record CFplus "CFplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CFplus,
    z=1);
   annotation (preferredView = "info");
  end CFplus;

  record CFBr3 "CFBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CFBr3);
   annotation (preferredView = "info");
  end CFBr3;

  record CFCL "CFCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCL);
   annotation (preferredView = "info");
  end CFCL;

  record CFCLBr2 "CFCLBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCLBr2);
   annotation (preferredView = "info");
  end CFCLBr2;

  record CFCL2 "CFCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCL2);
   annotation (preferredView = "info");
  end CFCL2;

  record CFCL2Br "CFCL2Br(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCL2Br);
   annotation (preferredView = "info");
  end CFCL2Br;

  record CFCL3 "CFCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCL3);
   annotation (preferredView = "info");
  end CFCL3;

  record CF2 "CF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2);
   annotation (preferredView = "info");
  end CF2;

  record CF2plus "CF2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2plus,
    z=1);
   annotation (preferredView = "info");
  end CF2plus;

  record CF2Br2 "CF2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2Br2);
   annotation (preferredView = "info");
  end CF2Br2;

  record CF2CL "CF2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2CL);
   annotation (preferredView = "info");
  end CF2CL;

  record CF2CLBr "CF2CLBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2CLBr);
   annotation (preferredView = "info");
  end CF2CLBr;

  record CF2CL2 "CF2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2CL2);
   annotation (preferredView = "info");
  end CF2CL2;

  record CF3 "CF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF3);
   annotation (preferredView = "info");
  end CF3;

  record CF3plus "CF3plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF3plus,
    z=1);
   annotation (preferredView = "info");
  end CF3plus;

  record CF3Br "CF3Br(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF3Br);
   annotation (preferredView = "info");
  end CF3Br;

  record CF3CL "CF3CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF3CL);
   annotation (preferredView = "info");
  end CF3CL;

  record CF4 "CF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CF4);
   annotation (preferredView = "info");
  end CF4;

  record CHplus "CHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHplus,
    z=1);
   annotation (preferredView = "info");
  end CHplus;

  record CHBr3 "CHBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHBr3);
   annotation (preferredView = "info");
  end CHBr3;

  record CHCL "CHCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCL);
   annotation (preferredView = "info");
  end CHCL;

  record CHCLBr2 "CHCLBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCLBr2);
   annotation (preferredView = "info");
  end CHCLBr2;

  record CHCL2 "CHCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCL2);
   annotation (preferredView = "info");
  end CHCL2;

  record CHCL2Br "CHCL2Br(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCL2Br);
   annotation (preferredView = "info");
  end CHCL2Br;

  record CHCL3 "CHCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCL3);
   annotation (preferredView = "info");
  end CHCL3;

  record CHF "CHF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF);
   annotation (preferredView = "info");
  end CHF;

  record CHFBr2 "CHFBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHFBr2);
   annotation (preferredView = "info");
  end CHFBr2;

  record CHFCL "CHFCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHFCL);
   annotation (preferredView = "info");
  end CHFCL;

  record CHFCLBr "CHFCLBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHFCLBr);
   annotation (preferredView = "info");
  end CHFCLBr;

  record CHFCL2 "CHFCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHFCL2);
   annotation (preferredView = "info");
  end CHFCL2;

  record CHF2 "CHF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF2);
   annotation (preferredView = "info");
  end CHF2;

  record CHF2Br "CHF2Br(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF2Br);
   annotation (preferredView = "info");
  end CHF2Br;

  record CHF2CL "CHF2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF2CL);
   annotation (preferredView = "info");
  end CHF2CL;

  record CHF3 "CHF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF3);
   annotation (preferredView = "info");
  end CHF3;

  record CHI3 "CHI3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CHI3);
   annotation (preferredView = "info");
  end CHI3;

  record CH2 "CH2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2);
   annotation (preferredView = "info");
  end CH2;

  record CH2Br2 "CH2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2Br2);
   annotation (preferredView = "info");
  end CH2Br2;

  record CH2CL "CH2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CL);
   annotation (preferredView = "info");
  end CH2CL;

  record CH2CLBr "CH2CLBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CLBr);
   annotation (preferredView = "info");
  end CH2CLBr;

  record CH2CL2 "CH2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CL2);
   annotation (preferredView = "info");
  end CH2CL2;

  record CH2F "CH2F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2F);
   annotation (preferredView = "info");
  end CH2F;

  record CH2FBr "CH2FBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2FBr);
   annotation (preferredView = "info");
  end CH2FBr;

  record CH2FCL "CH2FCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2FCL);
   annotation (preferredView = "info");
  end CH2FCL;

  record CH2F2 "CH2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2F2);
   annotation (preferredView = "info");
  end CH2F2;

  record CH2I2 "CH2I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2I2);
   annotation (preferredView = "info");
  end CH2I2;

  record CH3 "CH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3);
   annotation (preferredView = "info");
  end CH3;

  record CH3Br "CH3Br(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3Br);
   annotation (preferredView = "info");
  end CH3Br;

  record CH3CL "CH3CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3CL);
   annotation (preferredView = "info");
  end CH3CL;

  record CH3F "CH3F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3F);
   annotation (preferredView = "info");
  end CH3F;

  record CH3I "CH3I(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3I);
   annotation (preferredView = "info");
  end CH3I;

  record CH2OH "CH2OH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2OH);
   annotation (preferredView = "info");
  end CH2OH;

  record CH2OHplus "CH2OHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2OHplus,
    z=1);
   annotation (preferredView = "info");
  end CH2OHplus;

  record CH3O "CH3O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3O);
   annotation (preferredView = "info");
  end CH3O;

  record CH4 "CH4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH4);
   annotation (preferredView = "info");
  end CH4;

  record CH3OH "CH3OH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3OH);
   annotation (preferredView = "info");
  end CH3OH;

  record CH3OOH "CH3OOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3OOH);
   annotation (preferredView = "info");
  end CH3OOH;

  record CI "CI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CI);
   annotation (preferredView = "info");
  end CI;

  record CI2 "CI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CI2);
   annotation (preferredView = "info");
  end CI2;

  record CI3 "CI3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CI3);
   annotation (preferredView = "info");
  end CI3;

  record CI4 "CI4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CI4);
   annotation (preferredView = "info");
  end CI4;

  record CN "CN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CN);
   annotation (preferredView = "info");
  end CN;

  record CNplus "CNplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CNplus,
    z=1);
   annotation (preferredView = "info");
  end CNplus;

  record CNminus "CNminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CNminus,
    z=-1);
   annotation (preferredView = "info");
  end CNminus;

  record CNN "CNN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CNN);
   annotation (preferredView = "info");
  end CNN;

  record CO "CO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CO);
   annotation (preferredView = "info");
  end CO;

  record COplus "COplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.COplus,
    z=1);
   annotation (preferredView = "info");
  end COplus;

  record COCL "COCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.COCL);
   annotation (preferredView = "info");
  end COCL;

  record COCL2 "COCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.COCL2);
   annotation (preferredView = "info");
  end COCL2;

  record COFCL "COFCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.COFCL);
   annotation (preferredView = "info");
  end COFCL;

  record COF2 "COF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.COF2);
   annotation (preferredView = "info");
  end COF2;

  record COHCL "COHCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.COHCL);
   annotation (preferredView = "info");
  end COHCL;

  record COHF "COHF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.COHF);
   annotation (preferredView = "info");
  end COHF;

  record COS "COS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.COS);
   annotation (preferredView = "info");
  end COS;

  record CO2 "CO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CO2);
   annotation (preferredView = "info");
  end CO2;

  record CO2plus "CO2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CO2plus,
    z=1);
   annotation (preferredView = "info");
  end CO2plus;

  record COOH "COOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.COOH);
   annotation (preferredView = "info");
  end COOH;

  record CP "CP(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CP);
   annotation (preferredView = "info");
  end CP;

  record CS "CS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CS);
   annotation (preferredView = "info");
  end CS;

  record CS2 "CS2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CS2);
   annotation (preferredView = "info");
  end CS2;

  record C2 "C2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2);
   annotation (preferredView = "info");
  end C2;

  record C2plus "C2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2plus,
    z=1);
   annotation (preferredView = "info");
  end C2plus;

  record C2minus "C2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2minus,
    z=-1);
   annotation (preferredView = "info");
  end C2minus;

  record C2CL "C2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL);
   annotation (preferredView = "info");
  end C2CL;

  record C2CL2 "C2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL2);
   annotation (preferredView = "info");
  end C2CL2;

  record C2CL3 "C2CL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL3);
   annotation (preferredView = "info");
  end C2CL3;

  record C2CL4 "C2CL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL4);
   annotation (preferredView = "info");
  end C2CL4;

  record C2CL6 "C2CL6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL6);
   annotation (preferredView = "info");
  end C2CL6;

  record C2F "C2F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F);
   annotation (preferredView = "info");
  end C2F;

  record C2FCL "C2FCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2FCL);
   annotation (preferredView = "info");
  end C2FCL;

  record C2FCL3 "C2FCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2FCL3);
   annotation (preferredView = "info");
  end C2FCL3;

  record C2F2 "C2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F2);
   annotation (preferredView = "info");
  end C2F2;

  record C2F2CL2 "C2F2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F2CL2);
   annotation (preferredView = "info");
  end C2F2CL2;

  record C2F3 "C2F3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F3);
   annotation (preferredView = "info");
  end C2F3;

  record C2F3CL "C2F3CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F3CL);
   annotation (preferredView = "info");
  end C2F3CL;

  record C2F4 "C2F4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F4);
   annotation (preferredView = "info");
  end C2F4;

  record C2F6 "C2F6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F6);
   annotation (preferredView = "info");
  end C2F6;

  record C2H "C2H(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H);
   annotation (preferredView = "info");
  end C2H;

  record C2HCL "C2HCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HCL);
   annotation (preferredView = "info");
  end C2HCL;

  record C2HCL3 "C2HCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HCL3);
   annotation (preferredView = "info");
  end C2HCL3;

  record C2HF "C2HF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HF);
   annotation (preferredView = "info");
  end C2HF;

  record C2HFCL2 "C2HFCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HFCL2);
   annotation (preferredView = "info");
  end C2HFCL2;

  record C2HF2CL "C2HF2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HF2CL);
   annotation (preferredView = "info");
  end C2HF2CL;

  record C2HF3 "C2HF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HF3);
   annotation (preferredView = "info");
  end C2HF3;

  record C2H2_vinylidene "C2H2_vinylidene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H2_vinylidene);
   annotation (preferredView = "info");
  end C2H2_vinylidene;

  record C2H2CL2 "C2H2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H2CL2);
   annotation (preferredView = "info");
  end C2H2CL2;

  record C2H2FCL "C2H2FCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H2FCL);
   annotation (preferredView = "info");
  end C2H2FCL;

  record C2H2F2 "C2H2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H2F2);
   annotation (preferredView = "info");
  end C2H2F2;

  record CH2CO_ketene "CH2CO_ketene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CO_ketene);
   annotation (preferredView = "info");
  end CH2CO_ketene;

  record O_CH_2O "O_CH_2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.O_CH_2O);
   annotation (preferredView = "info");
  end O_CH_2O;

  record HO_CO_2OH "HO_CO_2OH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HO_CO_2OH);
   annotation (preferredView = "info");
  end HO_CO_2OH;

  record C2H3_vinyl "C2H3_vinyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H3_vinyl);
   annotation (preferredView = "info");
  end C2H3_vinyl;

  record CH2BrminusCOOH "CH2BrminusCOOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2BrminusCOOH);
   annotation (preferredView = "info");
  end CH2BrminusCOOH;

  record C2H3CL "C2H3CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H3CL);
   annotation (preferredView = "info");
  end C2H3CL;

  record CH2CLminusCOOH "CH2CLminusCOOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CLminusCOOH);
   annotation (preferredView = "info");
  end CH2CLminusCOOH;

  record C2H3F "C2H3F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H3F);
   annotation (preferredView = "info");
  end C2H3F;

  record CH3CN "CH3CN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3CN);
   annotation (preferredView = "info");
  end CH3CN;

  record CH3CO_acetyl "CH3CO_acetyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3CO_acetyl);
   annotation (preferredView = "info");
  end CH3CO_acetyl;

  record C2H4 "C2H4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H4);
   annotation (preferredView = "info");
  end C2H4;

  record C2H4O_ethylen_o "C2H4O_ethylen_o(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H4O_ethylen_o);
   annotation (preferredView = "info");
  end C2H4O_ethylen_o;

  record CH3CHO_ethanal "CH3CHO_ethanal(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3CHO_ethanal);
   annotation (preferredView = "info");
  end CH3CHO_ethanal;

  record CH3COOH "CH3COOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3COOH);
   annotation (preferredView = "info");
  end CH3COOH;

  record OHCH2COOH "OHCH2COOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.OHCH2COOH);
   annotation (preferredView = "info");
  end OHCH2COOH;

  record C2H5 "C2H5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H5);
   annotation (preferredView = "info");
  end C2H5;

  record C2H5Br "C2H5Br(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H5Br);
   annotation (preferredView = "info");
  end C2H5Br;

  record C2H6 "C2H6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H6);
   annotation (preferredView = "info");
  end C2H6;

  record CH3N2CH3 "CH3N2CH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3N2CH3);
   annotation (preferredView = "info");
  end CH3N2CH3;

  record C2H5OH "C2H5OH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H5OH);
   annotation (preferredView = "info");
  end C2H5OH;

  record CH3OCH3 "CH3OCH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3OCH3);
   annotation (preferredView = "info");
  end CH3OCH3;

  record CH3O2CH3 "CH3O2CH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3O2CH3);
   annotation (preferredView = "info");
  end CH3O2CH3;

  record CCN "CCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CCN);
   annotation (preferredView = "info");
  end CCN;

  record CNC "CNC(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CNC);
   annotation (preferredView = "info");
  end CNC;

  record OCCN "OCCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.OCCN);
   annotation (preferredView = "info");
  end OCCN;

  record C2N2 "C2N2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2N2);
   annotation (preferredView = "info");
  end C2N2;

  record C2O "C2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C2O);
   annotation (preferredView = "info");
  end C2O;

  record C3 "C3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3);
   annotation (preferredView = "info");
  end C3;

  record C3H3_1_propynl "C3H3_1_propynl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H3_1_propynl);
   annotation (preferredView = "info");
  end C3H3_1_propynl;

  record C3H3_2_propynl "C3H3_2_propynl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H3_2_propynl);
   annotation (preferredView = "info");
  end C3H3_2_propynl;

  record C3H4_allene "C3H4_allene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H4_allene);
   annotation (preferredView = "info");
  end C3H4_allene;

  record C3H4_propyne "C3H4_propyne(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H4_propyne);
   annotation (preferredView = "info");
  end C3H4_propyne;

  record C3H4_cyclo "C3H4_cyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H4_cyclo);
   annotation (preferredView = "info");
  end C3H4_cyclo;

  record C3H5_allyl "C3H5_allyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H5_allyl);
   annotation (preferredView = "info");
  end C3H5_allyl;

  record C3H6_propylene "C3H6_propylene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6_propylene);
   annotation (preferredView = "info");
  end C3H6_propylene;

  record C3H6_cyclo "C3H6_cyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6_cyclo);
   annotation (preferredView = "info");
  end C3H6_cyclo;

  record C3H6O_propylox "C3H6O_propylox(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6O_propylox);
   annotation (preferredView = "info");
  end C3H6O_propylox;

  record C3H6O_acetone "C3H6O_acetone(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6O_acetone);
   annotation (preferredView = "info");
  end C3H6O_acetone;

  record C3H6O_propanal "C3H6O_propanal(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6O_propanal);
   annotation (preferredView = "info");
  end C3H6O_propanal;

  record C3H7_n_propyl "C3H7_n_propyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H7_n_propyl);
   annotation (preferredView = "info");
  end C3H7_n_propyl;

  record C3H7_i_propyl "C3H7_i_propyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H7_i_propyl);
   annotation (preferredView = "info");
  end C3H7_i_propyl;

  record C3H8 "C3H8(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H8);
   annotation (preferredView = "info");
  end C3H8;

  record C3H8O_1propanol "C3H8O_1propanol(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H8O_1propanol);
   annotation (preferredView = "info");
  end C3H8O_1propanol;

  record C3H8O_2propanol "C3H8O_2propanol(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H8O_2propanol);
   annotation (preferredView = "info");
  end C3H8O_2propanol;

  record CNCOCN "CNCOCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CNCOCN);
   annotation (preferredView = "info");
  end CNCOCN;

  record C3O2 "C3O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C3O2);
   annotation (preferredView = "info");
  end C3O2;

  record C4 "C4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4);
   annotation (preferredView = "info");
  end C4;

  record C4H2_butadiyne "C4H2_butadiyne(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H2_butadiyne);
   annotation (preferredView = "info");
  end C4H2_butadiyne;

  record C4H4_1_3minuscyclo "C4H4_1_3minuscyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H4_1_3minuscyclo);
   annotation (preferredView = "info");
  end C4H4_1_3minuscyclo;

  record C4H6_butadiene "C4H6_butadiene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H6_butadiene);
   annotation (preferredView = "info");
  end C4H6_butadiene;

  record C4H6_1butyne "C4H6_1butyne(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H6_1butyne);
   annotation (preferredView = "info");
  end C4H6_1butyne;

  record C4H6_2butyne "C4H6_2butyne(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H6_2butyne);
   annotation (preferredView = "info");
  end C4H6_2butyne;

  record C4H6_cyclo "C4H6_cyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H6_cyclo);
   annotation (preferredView = "info");
  end C4H6_cyclo;

  record C4H8_1_butene "C4H8_1_butene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H8_1_butene);
   annotation (preferredView = "info");
  end C4H8_1_butene;

  record C4H8_cis2_buten "C4H8_cis2_buten(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H8_cis2_buten);
   annotation (preferredView = "info");
  end C4H8_cis2_buten;

  record C4H8_isobutene "C4H8_isobutene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H8_isobutene);
   annotation (preferredView = "info");
  end C4H8_isobutene;

  record C4H8_cyclo "C4H8_cyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H8_cyclo);
   annotation (preferredView = "info");
  end C4H8_cyclo;

  record C4H9_n_butyl "C4H9_n_butyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H9_n_butyl);
   annotation (preferredView = "info");
  end C4H9_n_butyl;

  record C4H9_i_butyl "C4H9_i_butyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H9_i_butyl);
   annotation (preferredView = "info");
  end C4H9_i_butyl;

  record C4H9_s_butyl "C4H9_s_butyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H9_s_butyl);
   annotation (preferredView = "info");
  end C4H9_s_butyl;

  record C4H9_t_butyl "C4H9_t_butyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H9_t_butyl);
   annotation (preferredView = "info");
  end C4H9_t_butyl;

  record C4H10_n_butane "C4H10_n_butane(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H10_n_butane);
   annotation (preferredView = "info");
  end C4H10_n_butane;

  record C4H10_isobutane "C4H10_isobutane(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H10_isobutane);
   annotation (preferredView = "info");
  end C4H10_isobutane;

  record C4N2 "C4N2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C4N2);
   annotation (preferredView = "info");
  end C4N2;

  record C5 "C5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C5);
   annotation (preferredView = "info");
  end C5;

  record C5H6_1_3cyclo "C5H6_1_3cyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H6_1_3cyclo);
   annotation (preferredView = "info");
  end C5H6_1_3cyclo;

  record C5H8_cyclo "C5H8_cyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H8_cyclo);
   annotation (preferredView = "info");
  end C5H8_cyclo;

  record C5H10_1_pentene "C5H10_1_pentene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H10_1_pentene);
   annotation (preferredView = "info");
  end C5H10_1_pentene;

  record C5H10_cyclo "C5H10_cyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H10_cyclo);
   annotation (preferredView = "info");
  end C5H10_cyclo;

  record C5H11_pentyl "C5H11_pentyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H11_pentyl);
   annotation (preferredView = "info");
  end C5H11_pentyl;

  record C5H11_t_pentyl "C5H11_t_pentyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H11_t_pentyl);
   annotation (preferredView = "info");
  end C5H11_t_pentyl;

  record C5H12_n_pentane "C5H12_n_pentane(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H12_n_pentane);
   annotation (preferredView = "info");
  end C5H12_n_pentane;

  record C5H12_i_pentane "C5H12_i_pentane(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H12_i_pentane);
   annotation (preferredView = "info");
  end C5H12_i_pentane;

  record CH3C_CH3_2CH3 "CH3C_CH3_2CH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3C_CH3_2CH3);
   annotation (preferredView = "info");
  end CH3C_CH3_2CH3;

  record C6D5_phenyl "C6D5_phenyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6D5_phenyl);
   annotation (preferredView = "info");
  end C6D5_phenyl;

  record C6D6 "C6D6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6D6);
   annotation (preferredView = "info");
  end C6D6;

  record C6H2 "C6H2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H2);
   annotation (preferredView = "info");
  end C6H2;

  record C6H5_phenyl "C6H5_phenyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H5_phenyl);
   annotation (preferredView = "info");
  end C6H5_phenyl;

  record C6H5O_phenoxy "C6H5O_phenoxy(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H5O_phenoxy);
   annotation (preferredView = "info");
  end C6H5O_phenoxy;

  record C6H6 "C6H6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H6);
   annotation (preferredView = "info");
  end C6H6;

  record C6H5OH_phenol "C6H5OH_phenol(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H5OH_phenol);
   annotation (preferredView = "info");
  end C6H5OH_phenol;

  record C6H10_cyclo "C6H10_cyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H10_cyclo);
   annotation (preferredView = "info");
  end C6H10_cyclo;

  record C6H12_1_hexene "C6H12_1_hexene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H12_1_hexene);
   annotation (preferredView = "info");
  end C6H12_1_hexene;

  record C6H12_cyclo "C6H12_cyclo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H12_cyclo);
   annotation (preferredView = "info");
  end C6H12_cyclo;

  record C6H13_n_hexyl "C6H13_n_hexyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H13_n_hexyl);
   annotation (preferredView = "info");
  end C6H13_n_hexyl;

  record C6H14_n_hexane "C6H14_n_hexane(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H14_n_hexane);
   annotation (preferredView = "info");
  end C6H14_n_hexane;

  record C7H7_benzyl "C7H7_benzyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H7_benzyl);
   annotation (preferredView = "info");
  end C7H7_benzyl;

  record C7H8 "C7H8(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H8);
   annotation (preferredView = "info");
  end C7H8;

  record C7H8O_cresol_mx "C7H8O_cresol_mx(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H8O_cresol_mx);
   annotation (preferredView = "info");
  end C7H8O_cresol_mx;

  record C7H14_1_heptene "C7H14_1_heptene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H14_1_heptene);
   annotation (preferredView = "info");
  end C7H14_1_heptene;

  record C7H15_n_heptyl "C7H15_n_heptyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H15_n_heptyl);
   annotation (preferredView = "info");
  end C7H15_n_heptyl;

  record C7H16_n_heptane "C7H16_n_heptane(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H16_n_heptane);
   annotation (preferredView = "info");
  end C7H16_n_heptane;

  record C7H16_2_methylh "C7H16_2_methylh(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H16_2_methylh);
   annotation (preferredView = "info");
  end C7H16_2_methylh;

  record C8H8_styrene "C8H8_styrene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H8_styrene);
   annotation (preferredView = "info");
  end C8H8_styrene;

  record C8H10_ethylbenz "C8H10_ethylbenz(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H10_ethylbenz);
   annotation (preferredView = "info");
  end C8H10_ethylbenz;

  record C8H16_1_octene "C8H16_1_octene(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H16_1_octene);
   annotation (preferredView = "info");
  end C8H16_1_octene;

  record C8H17_n_octyl "C8H17_n_octyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H17_n_octyl);
   annotation (preferredView = "info");
  end C8H17_n_octyl;

  record C8H18_n_octane "C8H18_n_octane(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H18_n_octane);
   annotation (preferredView = "info");
  end C8H18_n_octane;

  record C8H18_isooctane "C8H18_isooctane(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H18_isooctane);
   annotation (preferredView = "info");
  end C8H18_isooctane;

  record C9H19_n_nonyl "C9H19_n_nonyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C9H19_n_nonyl);
   annotation (preferredView = "info");
  end C9H19_n_nonyl;

  record C10H8_naphthale "C10H8_naphthale(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C10H8_naphthale);
   annotation (preferredView = "info");
  end C10H8_naphthale;

  record C10H21_n_decyl "C10H21_n_decyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C10H21_n_decyl);
   annotation (preferredView = "info");
  end C10H21_n_decyl;

  record C12H9_o_bipheny "C12H9_o_bipheny(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C12H9_o_bipheny);
   annotation (preferredView = "info");
  end C12H9_o_bipheny;

  record C12H10_biphenyl "C12H10_biphenyl(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.C12H10_biphenyl);
   annotation (preferredView = "info");
  end C12H10_biphenyl;

  record Ca "Ca(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ca);
   annotation (preferredView = "info");
  end Ca;

  record Caplus "Caplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Caplus,
    z=1);
   annotation (preferredView = "info");
  end Caplus;

  record CaBr "CaBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaBr);
   annotation (preferredView = "info");
  end CaBr;

  record CaBr2 "CaBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaBr2);
   annotation (preferredView = "info");
  end CaBr2;

  record CaCL "CaCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaCL);
   annotation (preferredView = "info");
  end CaCL;

  record CaCLplus "CaCLplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaCLplus,
    z=1);
   annotation (preferredView = "info");
  end CaCLplus;

  record CaCL2 "CaCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaCL2);
   annotation (preferredView = "info");
  end CaCL2;

  record CaF "CaF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaF);
   annotation (preferredView = "info");
  end CaF;

  record CaFplus "CaFplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaFplus,
    z=1);
   annotation (preferredView = "info");
  end CaFplus;

  record CaF2 "CaF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaF2);
   annotation (preferredView = "info");
  end CaF2;

  record CaH "CaH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaH);
   annotation (preferredView = "info");
  end CaH;

  record CaI "CaI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaI);
   annotation (preferredView = "info");
  end CaI;

  record CaI2 "CaI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaI2);
   annotation (preferredView = "info");
  end CaI2;

  record CaO "CaO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaO);
   annotation (preferredView = "info");
  end CaO;

  record CaOplus "CaOplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaOplus,
    z=1);
   annotation (preferredView = "info");
  end CaOplus;

  record CaOH "CaOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaOH);
   annotation (preferredView = "info");
  end CaOH;

  record CaOHplus "CaOHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaOHplus,
    z=1);
   annotation (preferredView = "info");
  end CaOHplus;

  record Ca_OH_2 "Ca_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ca_OH_2);
   annotation (preferredView = "info");
  end Ca_OH_2;

  record CaS "CaS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CaS);
   annotation (preferredView = "info");
  end CaS;

  record Ca2 "Ca2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ca2);
   annotation (preferredView = "info");
  end Ca2;

  record Cd "Cd(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cd);
   annotation (preferredView = "info");
  end Cd;

  record Cdplus "Cdplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cdplus,
    z=1);
   annotation (preferredView = "info");
  end Cdplus;

  record CL "CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CL);
   annotation (preferredView = "info");
  end CL;

  record CLplus "CLplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CLplus,
    z=1);
   annotation (preferredView = "info");
  end CLplus;

  record CLminus "CLminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CLminus,
    z=-1);
   annotation (preferredView = "info");
  end CLminus;

  record CLCN "CLCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CLCN);
   annotation (preferredView = "info");
  end CLCN;

  record CLF "CLF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CLF);
   annotation (preferredView = "info");
  end CLF;

  record CLF3 "CLF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CLF3);
   annotation (preferredView = "info");
  end CLF3;

  record CLF5 "CLF5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CLF5);
   annotation (preferredView = "info");
  end CLF5;

  record CLO "CLO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CLO);
   annotation (preferredView = "info");
  end CLO;

  record CLO2 "CLO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CLO2);
   annotation (preferredView = "info");
  end CLO2;

  record CL2 "CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CL2);
   annotation (preferredView = "info");
  end CL2;

  record CL2O "CL2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CL2O);
   annotation (preferredView = "info");
  end CL2O;

  record Co "Co(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Co);
   annotation (preferredView = "info");
  end Co;

  record Coplus "Coplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Coplus,
    z=1);
   annotation (preferredView = "info");
  end Coplus;

  record Cominus "Cominus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cominus,
    z=-1);
   annotation (preferredView = "info");
  end Cominus;

  record Cr "Cr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cr);
   annotation (preferredView = "info");
  end Cr;

  record Crplus "Crplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Crplus,
    z=1);
   annotation (preferredView = "info");
  end Crplus;

  record Crminus "Crminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Crminus,
    z=-1);
   annotation (preferredView = "info");
  end Crminus;

  record CrN "CrN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CrN);
   annotation (preferredView = "info");
  end CrN;

  record CrO "CrO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CrO);
   annotation (preferredView = "info");
  end CrO;

  record CrO2 "CrO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CrO2);
   annotation (preferredView = "info");
  end CrO2;

  record CrO3 "CrO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CrO3);
   annotation (preferredView = "info");
  end CrO3;

  record CrO3minus "CrO3minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CrO3minus,
    z=-1);
   annotation (preferredView = "info");
  end CrO3minus;

  record Cs "Cs(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs);
   annotation (preferredView = "info");
  end Cs;

  record Csplus "Csplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Csplus,
    z=1);
   annotation (preferredView = "info");
  end Csplus;

  record Csminus "Csminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Csminus,
    z=-1);
   annotation (preferredView = "info");
  end Csminus;

  record CsBO2 "CsBO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsBO2);
   annotation (preferredView = "info");
  end CsBO2;

  record CsBr "CsBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsBr);
   annotation (preferredView = "info");
  end CsBr;

  record CsCL "CsCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsCL);
   annotation (preferredView = "info");
  end CsCL;

  record CsF "CsF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsF);
   annotation (preferredView = "info");
  end CsF;

  record CsH "CsH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsH);
   annotation (preferredView = "info");
  end CsH;

  record CsI "CsI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsI);
   annotation (preferredView = "info");
  end CsI;

  record CsLi "CsLi(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsLi);
   annotation (preferredView = "info");
  end CsLi;

  record CsNO2 "CsNO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsNO2);
   annotation (preferredView = "info");
  end CsNO2;

  record CsNO3 "CsNO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsNO3);
   annotation (preferredView = "info");
  end CsNO3;

  record CsNa "CsNa(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsNa);
   annotation (preferredView = "info");
  end CsNa;

  record CsO "CsO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsO);
   annotation (preferredView = "info");
  end CsO;

  record CsOH "CsOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsOH);
   annotation (preferredView = "info");
  end CsOH;

  record CsRb "CsRb(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CsRb);
   annotation (preferredView = "info");
  end CsRb;

  record Cs2 "Cs2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2);
   annotation (preferredView = "info");
  end Cs2;

  record Cs2Br2 "Cs2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2Br2);
   annotation (preferredView = "info");
  end Cs2Br2;

  record Cs2CO3 "Cs2CO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2CO3);
   annotation (preferredView = "info");
  end Cs2CO3;

  record Cs2CL2 "Cs2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2CL2);
   annotation (preferredView = "info");
  end Cs2CL2;

  record Cs2F2 "Cs2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2F2);
   annotation (preferredView = "info");
  end Cs2F2;

  record Cs2I2 "Cs2I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2I2);
   annotation (preferredView = "info");
  end Cs2I2;

  record Cs2O "Cs2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2O);
   annotation (preferredView = "info");
  end Cs2O;

  record Cs2Oplus "Cs2Oplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2Oplus,
    z=1);
   annotation (preferredView = "info");
  end Cs2Oplus;

  record Cs2O2 "Cs2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2O2);
   annotation (preferredView = "info");
  end Cs2O2;

  record Cs2O2H2 "Cs2O2H2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2O2H2);
   annotation (preferredView = "info");
  end Cs2O2H2;

  record Cs2SO4 "Cs2SO4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2SO4);
   annotation (preferredView = "info");
  end Cs2SO4;

  record Cu "Cu(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cu);
   annotation (preferredView = "info");
  end Cu;

  record Cuplus "Cuplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cuplus,
    z=1);
   annotation (preferredView = "info");
  end Cuplus;

  record Cuminus "Cuminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cuminus,
    z=-1);
   annotation (preferredView = "info");
  end Cuminus;

  record CuCL "CuCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CuCL);
   annotation (preferredView = "info");
  end CuCL;

  record CuF "CuF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CuF);
   annotation (preferredView = "info");
  end CuF;

  record CuF2 "CuF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CuF2);
   annotation (preferredView = "info");
  end CuF2;

  record CuO "CuO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.CuO);
   annotation (preferredView = "info");
  end CuO;

  record Cu2 "Cu2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cu2);
   annotation (preferredView = "info");
  end Cu2;

  record Cu3CL3 "Cu3CL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Cu3CL3);
   annotation (preferredView = "info");
  end Cu3CL3;

  record D "D(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.D);
   annotation (preferredView = "info");
  end D;

  record Dplus "Dplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Dplus,
    z=1);
   annotation (preferredView = "info");
  end Dplus;

  record Dminus "Dminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Dminus,
    z=-1);
   annotation (preferredView = "info");
  end Dminus;

  record DBr "DBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.DBr);
   annotation (preferredView = "info");
  end DBr;

  record DCL "DCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.DCL);
   annotation (preferredView = "info");
  end DCL;

  record DF "DF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.DF);
   annotation (preferredView = "info");
  end DF;

  record DOCL "DOCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.DOCL);
   annotation (preferredView = "info");
  end DOCL;

  record DO2 "DO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.DO2);
   annotation (preferredView = "info");
  end DO2;

  record DO2minus "DO2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.DO2minus,
    z=-1);
   annotation (preferredView = "info");
  end DO2minus;

  record D2 "D2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.D2);
   annotation (preferredView = "info");
  end D2;

  record D2plus "D2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.D2plus,
    z=1);
   annotation (preferredView = "info");
  end D2plus;

  record D2minus "D2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.D2minus,
    z=-1);
   annotation (preferredView = "info");
  end D2minus;

  record D2O "D2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.D2O);
   annotation (preferredView = "info");
  end D2O;

  record D2O2 "D2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.D2O2);
   annotation (preferredView = "info");
  end D2O2;

  record D2S "D2S(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.D2S);
   annotation (preferredView = "info");
  end D2S;

  record eminus "eminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.eminus,
    z=-1);
   annotation (preferredView = "info");
  end eminus;

  record F "F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.F);
   annotation (preferredView = "info");
  end F;

  record Fplus "Fplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Fplus,
    z=1);
   annotation (preferredView = "info");
  end Fplus;

  record Fminus "Fminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Fminus,
    z=-1);
   annotation (preferredView = "info");
  end Fminus;

  record FCN "FCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FCN);
   annotation (preferredView = "info");
  end FCN;

  record FCO "FCO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FCO);
   annotation (preferredView = "info");
  end FCO;

  record FO "FO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FO);
   annotation (preferredView = "info");
  end FO;

  record FO2_FOO "FO2_FOO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FO2_FOO);
   annotation (preferredView = "info");
  end FO2_FOO;

  record FO2_OFO "FO2_OFO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FO2_OFO);
   annotation (preferredView = "info");
  end FO2_OFO;

  record F2 "F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.F2);
   annotation (preferredView = "info");
  end F2;

  record F2O "F2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.F2O);
   annotation (preferredView = "info");
  end F2O;

  record F2O2 "F2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.F2O2);
   annotation (preferredView = "info");
  end F2O2;

  record FS2F "FS2F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FS2F);
   annotation (preferredView = "info");
  end FS2F;

  record Fe "Fe(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe);
   annotation (preferredView = "info");
  end Fe;

  record Feplus "Feplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Feplus,
    z=1);
   annotation (preferredView = "info");
  end Feplus;

  record Fe_CO_5 "Fe_CO_5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe_CO_5);
   annotation (preferredView = "info");
  end Fe_CO_5;

  record FeCL "FeCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FeCL);
   annotation (preferredView = "info");
  end FeCL;

  record FeCL2 "FeCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FeCL2);
   annotation (preferredView = "info");
  end FeCL2;

  record FeCL3 "FeCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FeCL3);
   annotation (preferredView = "info");
  end FeCL3;

  record FeO "FeO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.FeO);
   annotation (preferredView = "info");
  end FeO;

  record Fe_OH_2 "Fe_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe_OH_2);
   annotation (preferredView = "info");
  end Fe_OH_2;

  record Fe2CL4 "Fe2CL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe2CL4);
   annotation (preferredView = "info");
  end Fe2CL4;

  record Fe2CL6 "Fe2CL6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe2CL6);
   annotation (preferredView = "info");
  end Fe2CL6;

  record Ga "Ga(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga);
   annotation (preferredView = "info");
  end Ga;

  record Gaplus "Gaplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Gaplus,
    z=1);
   annotation (preferredView = "info");
  end Gaplus;

  record GaBr "GaBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaBr);
   annotation (preferredView = "info");
  end GaBr;

  record GaBr2 "GaBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaBr2);
   annotation (preferredView = "info");
  end GaBr2;

  record GaBr3 "GaBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaBr3);
   annotation (preferredView = "info");
  end GaBr3;

  record GaCL "GaCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaCL);
   annotation (preferredView = "info");
  end GaCL;

  record GaCL2 "GaCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaCL2);
   annotation (preferredView = "info");
  end GaCL2;

  record GaCL3 "GaCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaCL3);
   annotation (preferredView = "info");
  end GaCL3;

  record GaF "GaF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaF);
   annotation (preferredView = "info");
  end GaF;

  record GaF2 "GaF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaF2);
   annotation (preferredView = "info");
  end GaF2;

  record GaF3 "GaF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaF3);
   annotation (preferredView = "info");
  end GaF3;

  record GaH "GaH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaH);
   annotation (preferredView = "info");
  end GaH;

  record GaI "GaI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaI);
   annotation (preferredView = "info");
  end GaI;

  record GaI2 "GaI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaI2);
   annotation (preferredView = "info");
  end GaI2;

  record GaI3 "GaI3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaI3);
   annotation (preferredView = "info");
  end GaI3;

  record GaO "GaO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaO);
   annotation (preferredView = "info");
  end GaO;

  record GaOH "GaOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GaOH);
   annotation (preferredView = "info");
  end GaOH;

  record Ga2Br2 "Ga2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2Br2);
   annotation (preferredView = "info");
  end Ga2Br2;

  record Ga2Br4 "Ga2Br4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2Br4);
   annotation (preferredView = "info");
  end Ga2Br4;

  record Ga2Br6 "Ga2Br6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2Br6);
   annotation (preferredView = "info");
  end Ga2Br6;

  record Ga2CL2 "Ga2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2CL2);
   annotation (preferredView = "info");
  end Ga2CL2;

  record Ga2CL4 "Ga2CL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2CL4);
   annotation (preferredView = "info");
  end Ga2CL4;

  record Ga2CL6 "Ga2CL6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2CL6);
   annotation (preferredView = "info");
  end Ga2CL6;

  record Ga2F2 "Ga2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2F2);
   annotation (preferredView = "info");
  end Ga2F2;

  record Ga2F4 "Ga2F4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2F4);
   annotation (preferredView = "info");
  end Ga2F4;

  record Ga2F6 "Ga2F6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2F6);
   annotation (preferredView = "info");
  end Ga2F6;

  record Ga2I2 "Ga2I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2I2);
   annotation (preferredView = "info");
  end Ga2I2;

  record Ga2I4 "Ga2I4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2I4);
   annotation (preferredView = "info");
  end Ga2I4;

  record Ga2I6 "Ga2I6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2I6);
   annotation (preferredView = "info");
  end Ga2I6;

  record Ga2O "Ga2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2O);
   annotation (preferredView = "info");
  end Ga2O;

  record Ge "Ge(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ge);
   annotation (preferredView = "info");
  end Ge;

  record Geplus "Geplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Geplus,
    z=1);
   annotation (preferredView = "info");
  end Geplus;

  record Geminus "Geminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Geminus,
    z=-1);
   annotation (preferredView = "info");
  end Geminus;

  record GeBr "GeBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeBr);
   annotation (preferredView = "info");
  end GeBr;

  record GeBr2 "GeBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeBr2);
   annotation (preferredView = "info");
  end GeBr2;

  record GeBr3 "GeBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeBr3);
   annotation (preferredView = "info");
  end GeBr3;

  record GeBr4 "GeBr4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeBr4);
   annotation (preferredView = "info");
  end GeBr4;

  record GeCL "GeCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeCL);
   annotation (preferredView = "info");
  end GeCL;

  record GeCL2 "GeCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeCL2);
   annotation (preferredView = "info");
  end GeCL2;

  record GeCL3 "GeCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeCL3);
   annotation (preferredView = "info");
  end GeCL3;

  record GeCL4 "GeCL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeCL4);
   annotation (preferredView = "info");
  end GeCL4;

  record GeF "GeF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeF);
   annotation (preferredView = "info");
  end GeF;

  record GeF2 "GeF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeF2);
   annotation (preferredView = "info");
  end GeF2;

  record GeF3 "GeF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeF3);
   annotation (preferredView = "info");
  end GeF3;

  record GeF4 "GeF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeF4);
   annotation (preferredView = "info");
  end GeF4;

  record GeH4 "GeH4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeH4);
   annotation (preferredView = "info");
  end GeH4;

  record GeI "GeI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeI);
   annotation (preferredView = "info");
  end GeI;

  record GeO "GeO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeO);
   annotation (preferredView = "info");
  end GeO;

  record GeO2 "GeO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeO2);
   annotation (preferredView = "info");
  end GeO2;

  record GeS "GeS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeS);
   annotation (preferredView = "info");
  end GeS;

  record GeS2 "GeS2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.GeS2);
   annotation (preferredView = "info");
  end GeS2;

  record Ge2 "Ge2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ge2);
   annotation (preferredView = "info");
  end Ge2;

  record H "H(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H);
   annotation (preferredView = "info");
  end H;

  record Hplus "Hplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Hplus,
    z=1);
   annotation (preferredView = "info");
  end Hplus;

  record Hminus "Hminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Hminus,
    z=-1);
   annotation (preferredView = "info");
  end Hminus;

  record HALO "HALO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HALO);
   annotation (preferredView = "info");
  end HALO;

  record HALO2 "HALO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HALO2);
   annotation (preferredView = "info");
  end HALO2;

  record HBO "HBO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HBO);
   annotation (preferredView = "info");
  end HBO;

  record HBOplus "HBOplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HBOplus,
    z=1);
   annotation (preferredView = "info");
  end HBOplus;

  record HBO2 "HBO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HBO2);
   annotation (preferredView = "info");
  end HBO2;

  record HBS "HBS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HBS);
   annotation (preferredView = "info");
  end HBS;

  record HBSplus "HBSplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HBSplus,
    z=1);
   annotation (preferredView = "info");
  end HBSplus;

  record HCN "HCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HCN);
   annotation (preferredView = "info");
  end HCN;

  record HCO "HCO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HCO);
   annotation (preferredView = "info");
  end HCO;

  record HCOplus "HCOplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HCOplus,
    z=1);
   annotation (preferredView = "info");
  end HCOplus;

  record HCCN "HCCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HCCN);
   annotation (preferredView = "info");
  end HCCN;

  record HCCO "HCCO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HCCO);
   annotation (preferredView = "info");
  end HCCO;

  record HCL "HCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HCL);
   annotation (preferredView = "info");
  end HCL;

  record HD "HD(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HD);
   annotation (preferredView = "info");
  end HD;

  record HDplus "HDplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HDplus,
    z=1);
   annotation (preferredView = "info");
  end HDplus;

  record HDO "HDO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HDO);
   annotation (preferredView = "info");
  end HDO;

  record HDO2 "HDO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HDO2);
   annotation (preferredView = "info");
  end HDO2;

  record HF "HF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HF);
   annotation (preferredView = "info");
  end HF;

  record HI "HI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HI);
   annotation (preferredView = "info");
  end HI;

  record HNC "HNC(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HNC);
   annotation (preferredView = "info");
  end HNC;

  record HNCO "HNCO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HNCO);
   annotation (preferredView = "info");
  end HNCO;

  record HNO "HNO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HNO);
   annotation (preferredView = "info");
  end HNO;

  record HNO2 "HNO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HNO2);
   annotation (preferredView = "info");
  end HNO2;

  record HNO3 "HNO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HNO3);
   annotation (preferredView = "info");
  end HNO3;

  record HOCL "HOCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HOCL);
   annotation (preferredView = "info");
  end HOCL;

  record HOF "HOF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HOF);
   annotation (preferredView = "info");
  end HOF;

  record HO2 "HO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HO2);
   annotation (preferredView = "info");
  end HO2;

  record HO2minus "HO2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HO2minus,
    z=-1);
   annotation (preferredView = "info");
  end HO2minus;

  record HPO "HPO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HPO);
   annotation (preferredView = "info");
  end HPO;

  record HSO3F "HSO3F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HSO3F);
   annotation (preferredView = "info");
  end HSO3F;

  record H2 "H2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2);
   annotation (preferredView = "info");
  end H2;

  record H2plus "H2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2plus,
    z=1);
   annotation (preferredView = "info");
  end H2plus;

  record H2minus "H2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2minus,
    z=-1);
   annotation (preferredView = "info");
  end H2minus;

  record HBOH "HBOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HBOH);
   annotation (preferredView = "info");
  end HBOH;

  record HCOOH "HCOOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HCOOH);
   annotation (preferredView = "info");
  end HCOOH;

  record H2F2 "H2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2F2);
   annotation (preferredView = "info");
  end H2F2;

  record H2O "H2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
   annotation (preferredView = "info");
  end H2O;

  record H2Oplus "H2Oplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2Oplus,
    z=1);
   annotation (preferredView = "info");
  end H2Oplus;

  record H2O2 "H2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2O2);
   annotation (preferredView = "info");
  end H2O2;

  record H2S "H2S(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2S);
   annotation (preferredView = "info");
  end H2S;

  record H2SO4 "H2SO4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2SO4);
   annotation (preferredView = "info");
  end H2SO4;

  record H2BOH "H2BOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H2BOH);
   annotation (preferredView = "info");
  end H2BOH;

  record HB_OH_2 "HB_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HB_OH_2);
   annotation (preferredView = "info");
  end HB_OH_2;

  record H3BO3 "H3BO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H3BO3);
   annotation (preferredView = "info");
  end H3BO3;

  record H3B3O3 "H3B3O3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H3B3O3);
   annotation (preferredView = "info");
  end H3B3O3;

  record H3B3O6 "H3B3O6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H3B3O6);
   annotation (preferredView = "info");
  end H3B3O6;

  record H3F3 "H3F3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H3F3);
   annotation (preferredView = "info");
  end H3F3;

  record H3Oplus "H3Oplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H3Oplus,
    z=1);
   annotation (preferredView = "info");
  end H3Oplus;

  record H4F4 "H4F4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H4F4);
   annotation (preferredView = "info");
  end H4F4;

  record H5F5 "H5F5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H5F5);
   annotation (preferredView = "info");
  end H5F5;

  record H6F6 "H6F6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H6F6);
   annotation (preferredView = "info");
  end H6F6;

  record H7F7 "H7F7(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.H7F7);
   annotation (preferredView = "info");
  end H7F7;

  record He "He(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.He);
   annotation (preferredView = "info");
  end He;

  record Heplus "Heplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Heplus,
    z=1);
   annotation (preferredView = "info");
  end Heplus;

  record Hg "Hg(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Hg);
   annotation (preferredView = "info");
  end Hg;

  record Hgplus "Hgplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Hgplus,
    z=1);
   annotation (preferredView = "info");
  end Hgplus;

  record HgBr2 "HgBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.HgBr2);
   annotation (preferredView = "info");
  end HgBr2;

  record I "I(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.I);
   annotation (preferredView = "info");
  end I;

  record Iplus "Iplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Iplus,
    z=1);
   annotation (preferredView = "info");
  end Iplus;

  record Iminus "Iminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Iminus,
    z=-1);
   annotation (preferredView = "info");
  end Iminus;

  record IF5 "IF5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.IF5);
   annotation (preferredView = "info");
  end IF5;

  record IF7 "IF7(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.IF7);
   annotation (preferredView = "info");
  end IF7;

  record I2 "I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.I2);
   annotation (preferredView = "info");
  end I2;

  record In "In(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In);
   annotation (preferredView = "info");
  end In;

  record Inplus "Inplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Inplus,
    z=1);
   annotation (preferredView = "info");
  end Inplus;

  record InBr "InBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InBr);
   annotation (preferredView = "info");
  end InBr;

  record InBr2 "InBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InBr2);
   annotation (preferredView = "info");
  end InBr2;

  record InBr3 "InBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InBr3);
   annotation (preferredView = "info");
  end InBr3;

  record InCL "InCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InCL);
   annotation (preferredView = "info");
  end InCL;

  record InCL2 "InCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InCL2);
   annotation (preferredView = "info");
  end InCL2;

  record InCL3 "InCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InCL3);
   annotation (preferredView = "info");
  end InCL3;

  record InF "InF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InF);
   annotation (preferredView = "info");
  end InF;

  record InF2 "InF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InF2);
   annotation (preferredView = "info");
  end InF2;

  record InF3 "InF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InF3);
   annotation (preferredView = "info");
  end InF3;

  record InH "InH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InH);
   annotation (preferredView = "info");
  end InH;

  record InI "InI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InI);
   annotation (preferredView = "info");
  end InI;

  record InI2 "InI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InI2);
   annotation (preferredView = "info");
  end InI2;

  record InI3 "InI3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InI3);
   annotation (preferredView = "info");
  end InI3;

  record InO "InO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InO);
   annotation (preferredView = "info");
  end InO;

  record InOH "InOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.InOH);
   annotation (preferredView = "info");
  end InOH;

  record In2Br2 "In2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2Br2);
   annotation (preferredView = "info");
  end In2Br2;

  record In2Br4 "In2Br4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2Br4);
   annotation (preferredView = "info");
  end In2Br4;

  record In2Br6 "In2Br6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2Br6);
   annotation (preferredView = "info");
  end In2Br6;

  record In2CL2 "In2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2CL2);
   annotation (preferredView = "info");
  end In2CL2;

  record In2CL4 "In2CL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2CL4);
   annotation (preferredView = "info");
  end In2CL4;

  record In2CL6 "In2CL6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2CL6);
   annotation (preferredView = "info");
  end In2CL6;

  record In2F2 "In2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2F2);
   annotation (preferredView = "info");
  end In2F2;

  record In2F4 "In2F4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2F4);
   annotation (preferredView = "info");
  end In2F4;

  record In2F6 "In2F6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2F6);
   annotation (preferredView = "info");
  end In2F6;

  record In2I2 "In2I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2I2);
   annotation (preferredView = "info");
  end In2I2;

  record In2I4 "In2I4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2I4);
   annotation (preferredView = "info");
  end In2I4;

  record In2I6 "In2I6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2I6);
   annotation (preferredView = "info");
  end In2I6;

  record In2O "In2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.In2O);
   annotation (preferredView = "info");
  end In2O;

  record K "K(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K);
   annotation (preferredView = "info");
  end K;

  record Kplus "Kplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Kplus,
    z=1);
   annotation (preferredView = "info");
  end Kplus;

  record Kminus "Kminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Kminus,
    z=-1);
   annotation (preferredView = "info");
  end Kminus;

  record KALF4 "KALF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KALF4);
   annotation (preferredView = "info");
  end KALF4;

  record KBO2 "KBO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KBO2);
   annotation (preferredView = "info");
  end KBO2;

  record KBr "KBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KBr);
   annotation (preferredView = "info");
  end KBr;

  record KCN "KCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KCN);
   annotation (preferredView = "info");
  end KCN;

  record KCL "KCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KCL);
   annotation (preferredView = "info");
  end KCL;

  record KF "KF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KF);
   annotation (preferredView = "info");
  end KF;

  record KH "KH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KH);
   annotation (preferredView = "info");
  end KH;

  record KI "KI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KI);
   annotation (preferredView = "info");
  end KI;

  record KLi "KLi(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KLi);
   annotation (preferredView = "info");
  end KLi;

  record KNO2 "KNO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KNO2);
   annotation (preferredView = "info");
  end KNO2;

  record KNO3 "KNO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KNO3);
   annotation (preferredView = "info");
  end KNO3;

  record KNa "KNa(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KNa);
   annotation (preferredView = "info");
  end KNa;

  record KO "KO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KO);
   annotation (preferredView = "info");
  end KO;

  record KOH "KOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.KOH);
   annotation (preferredView = "info");
  end KOH;

  record K2 "K2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2);
   annotation (preferredView = "info");
  end K2;

  record K2plus "K2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2plus,
    z=1);
   annotation (preferredView = "info");
  end K2plus;

  record K2Br2 "K2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2Br2);
   annotation (preferredView = "info");
  end K2Br2;

  record K2CO3 "K2CO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2CO3);
   annotation (preferredView = "info");
  end K2CO3;

  record K2C2N2 "K2C2N2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2C2N2);
   annotation (preferredView = "info");
  end K2C2N2;

  record K2CL2 "K2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2CL2);
   annotation (preferredView = "info");
  end K2CL2;

  record K2F2 "K2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2F2);
   annotation (preferredView = "info");
  end K2F2;

  record K2I2 "K2I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2I2);
   annotation (preferredView = "info");
  end K2I2;

  record K2O "K2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2O);
   annotation (preferredView = "info");
  end K2O;

  record K2Oplus "K2Oplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2Oplus,
    z=1);
   annotation (preferredView = "info");
  end K2Oplus;

  record K2O2 "K2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2O2);
   annotation (preferredView = "info");
  end K2O2;

  record K2O2H2 "K2O2H2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2O2H2);
   annotation (preferredView = "info");
  end K2O2H2;

  record K2SO4 "K2SO4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.K2SO4);
   annotation (preferredView = "info");
  end K2SO4;

  record Kr "Kr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Kr);
   annotation (preferredView = "info");
  end Kr;

  record Krplus "Krplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Krplus,
    z=1);
   annotation (preferredView = "info");
  end Krplus;

  record Li "Li(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li);
   annotation (preferredView = "info");
  end Li;

  record Liplus "Liplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Liplus,
    z=1);
   annotation (preferredView = "info");
  end Liplus;

  record Liminus "Liminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Liminus,
    z=-1);
   annotation (preferredView = "info");
  end Liminus;

  record LiALF4 "LiALF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiALF4);
   annotation (preferredView = "info");
  end LiALF4;

  record LiBO2 "LiBO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiBO2);
   annotation (preferredView = "info");
  end LiBO2;

  record LiBr "LiBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiBr);
   annotation (preferredView = "info");
  end LiBr;

  record LiCL "LiCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiCL);
   annotation (preferredView = "info");
  end LiCL;

  record LiF "LiF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiF);
   annotation (preferredView = "info");
  end LiF;

  record LiH "LiH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiH);
   annotation (preferredView = "info");
  end LiH;

  record LiI "LiI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiI);
   annotation (preferredView = "info");
  end LiI;

  record LiN "LiN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiN);
   annotation (preferredView = "info");
  end LiN;

  record LiNO2 "LiNO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiNO2);
   annotation (preferredView = "info");
  end LiNO2;

  record LiNO3 "LiNO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiNO3);
   annotation (preferredView = "info");
  end LiNO3;

  record LiO "LiO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiO);
   annotation (preferredView = "info");
  end LiO;

  record LiOF "LiOF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiOF);
   annotation (preferredView = "info");
  end LiOF;

  record LiOH "LiOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiOH);
   annotation (preferredView = "info");
  end LiOH;

  record LiON "LiON(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.LiON);
   annotation (preferredView = "info");
  end LiON;

  record Li2 "Li2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2);
   annotation (preferredView = "info");
  end Li2;

  record Li2plus "Li2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2plus,
    z=1);
   annotation (preferredView = "info");
  end Li2plus;

  record Li2Br2 "Li2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2Br2);
   annotation (preferredView = "info");
  end Li2Br2;

  record Li2F2 "Li2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2F2);
   annotation (preferredView = "info");
  end Li2F2;

  record Li2I2 "Li2I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2I2);
   annotation (preferredView = "info");
  end Li2I2;

  record Li2O "Li2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2O);
   annotation (preferredView = "info");
  end Li2O;

  record Li2Oplus "Li2Oplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2Oplus,
    z=1);
   annotation (preferredView = "info");
  end Li2Oplus;

  record Li2O2 "Li2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2O2);
   annotation (preferredView = "info");
  end Li2O2;

  record Li2O2H2 "Li2O2H2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2O2H2);
   annotation (preferredView = "info");
  end Li2O2H2;

  record Li2SO4 "Li2SO4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2SO4);
   annotation (preferredView = "info");
  end Li2SO4;

  record Li3plus "Li3plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3plus,
    z=1);
   annotation (preferredView = "info");
  end Li3plus;

  record Li3Br3 "Li3Br3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3Br3);
   annotation (preferredView = "info");
  end Li3Br3;

  record Li3CL3 "Li3CL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3CL3);
   annotation (preferredView = "info");
  end Li3CL3;

  record Li3F3 "Li3F3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3F3);
   annotation (preferredView = "info");
  end Li3F3;

  record Li3I3 "Li3I3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3I3);
   annotation (preferredView = "info");
  end Li3I3;

  record Mg "Mg(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mg);
   annotation (preferredView = "info");
  end Mg;

  record Mgplus "Mgplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mgplus,
    z=1);
   annotation (preferredView = "info");
  end Mgplus;

  record MgBr "MgBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgBr);
   annotation (preferredView = "info");
  end MgBr;

  record MgBr2 "MgBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgBr2);
   annotation (preferredView = "info");
  end MgBr2;

  record MgCL "MgCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgCL);
   annotation (preferredView = "info");
  end MgCL;

  record MgCLplus "MgCLplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgCLplus,
    z=1);
   annotation (preferredView = "info");
  end MgCLplus;

  record MgCL2 "MgCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgCL2);
   annotation (preferredView = "info");
  end MgCL2;

  record MgF "MgF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgF);
   annotation (preferredView = "info");
  end MgF;

  record MgFplus "MgFplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgFplus,
    z=1);
   annotation (preferredView = "info");
  end MgFplus;

  record MgF2 "MgF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgF2);
   annotation (preferredView = "info");
  end MgF2;

  record MgF2plus "MgF2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgF2plus,
    z=1);
   annotation (preferredView = "info");
  end MgF2plus;

  record MgH "MgH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgH);
   annotation (preferredView = "info");
  end MgH;

  record MgI "MgI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgI);
   annotation (preferredView = "info");
  end MgI;

  record MgI2 "MgI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgI2);
   annotation (preferredView = "info");
  end MgI2;

  record MgN "MgN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgN);
   annotation (preferredView = "info");
  end MgN;

  record MgO "MgO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgO);
   annotation (preferredView = "info");
  end MgO;

  record MgOH "MgOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgOH);
   annotation (preferredView = "info");
  end MgOH;

  record MgOHplus "MgOHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgOHplus,
    z=1);
   annotation (preferredView = "info");
  end MgOHplus;

  record Mg_OH_2 "Mg_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mg_OH_2);
   annotation (preferredView = "info");
  end Mg_OH_2;

  record MgS "MgS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MgS);
   annotation (preferredView = "info");
  end MgS;

  record Mg2 "Mg2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mg2);
   annotation (preferredView = "info");
  end Mg2;

  record Mg2F4 "Mg2F4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mg2F4);
   annotation (preferredView = "info");
  end Mg2F4;

  record Mn "Mn(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mn);
   annotation (preferredView = "info");
  end Mn;

  record Mnplus "Mnplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mnplus,
    z=1);
   annotation (preferredView = "info");
  end Mnplus;

  record Mo "Mo(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo);
   annotation (preferredView = "info");
  end Mo;

  record Moplus "Moplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Moplus,
    z=1);
   annotation (preferredView = "info");
  end Moplus;

  record Mominus "Mominus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mominus,
    z=-1);
   annotation (preferredView = "info");
  end Mominus;

  record MoO "MoO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MoO);
   annotation (preferredView = "info");
  end MoO;

  record MoO2 "MoO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MoO2);
   annotation (preferredView = "info");
  end MoO2;

  record MoO3 "MoO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MoO3);
   annotation (preferredView = "info");
  end MoO3;

  record MoO3minus "MoO3minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.MoO3minus,
    z=-1);
   annotation (preferredView = "info");
  end MoO3minus;

  record Mo2O6 "Mo2O6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo2O6);
   annotation (preferredView = "info");
  end Mo2O6;

  record Mo3O9 "Mo3O9(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo3O9);
   annotation (preferredView = "info");
  end Mo3O9;

  record Mo4O12 "Mo4O12(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo4O12);
   annotation (preferredView = "info");
  end Mo4O12;

  record Mo5O15 "Mo5O15(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo5O15);
   annotation (preferredView = "info");
  end Mo5O15;

  record N "N(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N);
   annotation (preferredView = "info");
  end N;

  record Nplus "Nplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Nplus,
    z=1);
   annotation (preferredView = "info");
  end Nplus;

  record Nminus "Nminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Nminus,
    z=-1);
   annotation (preferredView = "info");
  end Nminus;

  record NCO "NCO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NCO);
   annotation (preferredView = "info");
  end NCO;

  record ND "ND(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ND);
   annotation (preferredView = "info");
  end ND;

  record ND2 "ND2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ND2);
   annotation (preferredView = "info");
  end ND2;

  record ND3 "ND3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ND3);
   annotation (preferredView = "info");
  end ND3;

  record NF "NF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NF);
   annotation (preferredView = "info");
  end NF;

  record NF2 "NF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NF2);
   annotation (preferredView = "info");
  end NF2;

  record NF3 "NF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NF3);
   annotation (preferredView = "info");
  end NF3;

  record NH "NH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NH);
   annotation (preferredView = "info");
  end NH;

  record NHplus "NHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NHplus,
    z=1);
   annotation (preferredView = "info");
  end NHplus;

  record NHF "NHF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NHF);
   annotation (preferredView = "info");
  end NHF;

  record NHF2 "NHF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NHF2);
   annotation (preferredView = "info");
  end NHF2;

  record NH2 "NH2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NH2);
   annotation (preferredView = "info");
  end NH2;

  record NH2F "NH2F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NH2F);
   annotation (preferredView = "info");
  end NH2F;

  record NH3 "NH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NH3);
   annotation (preferredView = "info");
  end NH3;

  record NH2OH "NH2OH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NH2OH);
   annotation (preferredView = "info");
  end NH2OH;

  record NH4plus "NH4plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NH4plus,
    z=1);
   annotation (preferredView = "info");
  end NH4plus;

  record NO "NO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NO);
   annotation (preferredView = "info");
  end NO;

  record NOCL "NOCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NOCL);
   annotation (preferredView = "info");
  end NOCL;

  record NOF "NOF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NOF);
   annotation (preferredView = "info");
  end NOF;

  record NOF3 "NOF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NOF3);
   annotation (preferredView = "info");
  end NOF3;

  record NO2 "NO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NO2);
   annotation (preferredView = "info");
  end NO2;

  record NO2minus "NO2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NO2minus,
    z=-1);
   annotation (preferredView = "info");
  end NO2minus;

  record NO2CL "NO2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NO2CL);
   annotation (preferredView = "info");
  end NO2CL;

  record NO2F "NO2F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NO2F);
   annotation (preferredView = "info");
  end NO2F;

  record NO3 "NO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NO3);
   annotation (preferredView = "info");
  end NO3;

  record NO3minus "NO3minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NO3minus,
    z=-1);
   annotation (preferredView = "info");
  end NO3minus;

  record NO3F "NO3F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NO3F);
   annotation (preferredView = "info");
  end NO3F;

  record N2 "N2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2);
   annotation (preferredView = "info");
  end N2;

  record N2plus "N2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2plus,
    z=1);
   annotation (preferredView = "info");
  end N2plus;

  record N2minus "N2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2minus,
    z=-1);
   annotation (preferredView = "info");
  end N2minus;

  record NCN "NCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NCN);
   annotation (preferredView = "info");
  end NCN;

  record N2D2_cis "N2D2_cis(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2D2_cis);
   annotation (preferredView = "info");
  end N2D2_cis;

  record N2F2 "N2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2F2);
   annotation (preferredView = "info");
  end N2F2;

  record N2F4 "N2F4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2F4);
   annotation (preferredView = "info");
  end N2F4;

  record N2H2 "N2H2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2H2);
   annotation (preferredView = "info");
  end N2H2;

  record NH2NO2 "NH2NO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NH2NO2);
   annotation (preferredView = "info");
  end NH2NO2;

  record N2H4 "N2H4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2H4);
   annotation (preferredView = "info");
  end N2H4;

  record N2O "N2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2O);
   annotation (preferredView = "info");
  end N2O;

  record N2Oplus "N2Oplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2Oplus,
    z=1);
   annotation (preferredView = "info");
  end N2Oplus;

  record N2O3 "N2O3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2O3);
   annotation (preferredView = "info");
  end N2O3;

  record N2O4 "N2O4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2O4);
   annotation (preferredView = "info");
  end N2O4;

  record N2O5 "N2O5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N2O5);
   annotation (preferredView = "info");
  end N2O5;

  record N3 "N3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N3);
   annotation (preferredView = "info");
  end N3;

  record N3H "N3H(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.N3H);
   annotation (preferredView = "info");
  end N3H;

  record Na "Na(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na);
   annotation (preferredView = "info");
  end Na;

  record Naplus "Naplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Naplus,
    z=1);
   annotation (preferredView = "info");
  end Naplus;

  record Naminus "Naminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Naminus,
    z=-1);
   annotation (preferredView = "info");
  end Naminus;

  record NaALF4 "NaALF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaALF4);
   annotation (preferredView = "info");
  end NaALF4;

  record NaBO2 "NaBO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaBO2);
   annotation (preferredView = "info");
  end NaBO2;

  record NaBr "NaBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaBr);
   annotation (preferredView = "info");
  end NaBr;

  record NaCN "NaCN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaCN);
   annotation (preferredView = "info");
  end NaCN;

  record NaCL "NaCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaCL);
   annotation (preferredView = "info");
  end NaCL;

  record NaF "NaF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaF);
   annotation (preferredView = "info");
  end NaF;

  record NaH "NaH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaH);
   annotation (preferredView = "info");
  end NaH;

  record NaI "NaI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaI);
   annotation (preferredView = "info");
  end NaI;

  record NaLi "NaLi(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaLi);
   annotation (preferredView = "info");
  end NaLi;

  record NaNO2 "NaNO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaNO2);
   annotation (preferredView = "info");
  end NaNO2;

  record NaNO3 "NaNO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaNO3);
   annotation (preferredView = "info");
  end NaNO3;

  record NaO "NaO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaO);
   annotation (preferredView = "info");
  end NaO;

  record NaOH "NaOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaOH);
   annotation (preferredView = "info");
  end NaOH;

  record NaOHplus "NaOHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NaOHplus,
    z=1);
   annotation (preferredView = "info");
  end NaOHplus;

  record Na2 "Na2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2);
   annotation (preferredView = "info");
  end Na2;

  record Na2Br2 "Na2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2Br2);
   annotation (preferredView = "info");
  end Na2Br2;

  record Na2CL2 "Na2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2CL2);
   annotation (preferredView = "info");
  end Na2CL2;

  record Na2F2 "Na2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2F2);
   annotation (preferredView = "info");
  end Na2F2;

  record Na2I2 "Na2I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2I2);
   annotation (preferredView = "info");
  end Na2I2;

  record Na2O "Na2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2O);
   annotation (preferredView = "info");
  end Na2O;

  record Na2Oplus "Na2Oplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2Oplus,
    z=1);
   annotation (preferredView = "info");
  end Na2Oplus;

  record Na2O2 "Na2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2O2);
   annotation (preferredView = "info");
  end Na2O2;

  record Na2O2H2 "Na2O2H2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2O2H2);
   annotation (preferredView = "info");
  end Na2O2H2;

  record Na2SO4 "Na2SO4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2SO4);
   annotation (preferredView = "info");
  end Na2SO4;

  record Na3CL3 "Na3CL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na3CL3);
   annotation (preferredView = "info");
  end Na3CL3;

  record Na3F3 "Na3F3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Na3F3);
   annotation (preferredView = "info");
  end Na3F3;

  record Nb "Nb(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Nb);
   annotation (preferredView = "info");
  end Nb;

  record Nbplus "Nbplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Nbplus,
    z=1);
   annotation (preferredView = "info");
  end Nbplus;

  record Nbminus "Nbminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Nbminus,
    z=-1);
   annotation (preferredView = "info");
  end Nbminus;

  record NbCL5 "NbCL5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NbCL5);
   annotation (preferredView = "info");
  end NbCL5;

  record NbO "NbO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NbO);
   annotation (preferredView = "info");
  end NbO;

  record NbOCL3 "NbOCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NbOCL3);
   annotation (preferredView = "info");
  end NbOCL3;

  record NbO2 "NbO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NbO2);
   annotation (preferredView = "info");
  end NbO2;

  record Ne "Ne(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ne);
   annotation (preferredView = "info");
  end Ne;

  record Neplus "Neplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Neplus,
    z=1);
   annotation (preferredView = "info");
  end Neplus;

  record Ni "Ni(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ni);
   annotation (preferredView = "info");
  end Ni;

  record Niplus "Niplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Niplus,
    z=1);
   annotation (preferredView = "info");
  end Niplus;

  record Niminus "Niminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Niminus,
    z=-1);
   annotation (preferredView = "info");
  end Niminus;

  record NiCL "NiCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NiCL);
   annotation (preferredView = "info");
  end NiCL;

  record NiCL2 "NiCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NiCL2);
   annotation (preferredView = "info");
  end NiCL2;

  record NiO "NiO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NiO);
   annotation (preferredView = "info");
  end NiO;

  record NiS "NiS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.NiS);
   annotation (preferredView = "info");
  end NiS;

  record O "O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.O);
   annotation (preferredView = "info");
  end O;

  record Oplus "Oplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Oplus,
    z=1);
   annotation (preferredView = "info");
  end Oplus;

  record Ominus "Ominus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ominus,
    z=-1);
   annotation (preferredView = "info");
  end Ominus;

  record OD "OD(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.OD);
   annotation (preferredView = "info");
  end OD;

  record ODminus "ODminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ODminus,
    z=-1);
   annotation (preferredView = "info");
  end ODminus;

  record OH "OH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.OH);
   annotation (preferredView = "info");
  end OH;

  record OHplus "OHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.OHplus,
    z=1);
   annotation (preferredView = "info");
  end OHplus;

  record OHminus "OHminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.OHminus,
    z=-1);
   annotation (preferredView = "info");
  end OHminus;

  record O2 "O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.O2);
   annotation (preferredView = "info");
  end O2;

  record O2plus "O2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.O2plus,
    z=1);
   annotation (preferredView = "info");
  end O2plus;

  record O2minus "O2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.O2minus,
    z=-1);
   annotation (preferredView = "info");
  end O2minus;

  record O3 "O3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.O3);
   annotation (preferredView = "info");
  end O3;

  record P "P(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P);
   annotation (preferredView = "info");
  end P;

  record Pplus "Pplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Pplus,
    z=1);
   annotation (preferredView = "info");
  end Pplus;

  record Pminus "Pminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Pminus,
    z=-1);
   annotation (preferredView = "info");
  end Pminus;

  record PCL "PCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL);
   annotation (preferredView = "info");
  end PCL;

  record PCL2 "PCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL2);
   annotation (preferredView = "info");
  end PCL2;

  record PCL2minus "PCL2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL2minus,
    z=-1);
   annotation (preferredView = "info");
  end PCL2minus;

  record PCL3 "PCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL3);
   annotation (preferredView = "info");
  end PCL3;

  record PCL5 "PCL5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL5);
   annotation (preferredView = "info");
  end PCL5;

  record PF "PF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PF);
   annotation (preferredView = "info");
  end PF;

  record PFplus "PFplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PFplus,
    z=1);
   annotation (preferredView = "info");
  end PFplus;

  record PFminus "PFminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PFminus,
    z=-1);
   annotation (preferredView = "info");
  end PFminus;

  record PFCL "PFCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PFCL);
   annotation (preferredView = "info");
  end PFCL;

  record PFCLminus "PFCLminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PFCLminus,
    z=-1);
   annotation (preferredView = "info");
  end PFCLminus;

  record PFCL2 "PFCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PFCL2);
   annotation (preferredView = "info");
  end PFCL2;

  record PFCL4 "PFCL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PFCL4);
   annotation (preferredView = "info");
  end PFCL4;

  record PF2 "PF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PF2);
   annotation (preferredView = "info");
  end PF2;

  record PF2minus "PF2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PF2minus,
    z=-1);
   annotation (preferredView = "info");
  end PF2minus;

  record PF2CL "PF2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PF2CL);
   annotation (preferredView = "info");
  end PF2CL;

  record PF2CL3 "PF2CL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PF2CL3);
   annotation (preferredView = "info");
  end PF2CL3;

  record PF3 "PF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PF3);
   annotation (preferredView = "info");
  end PF3;

  record PF3CL2 "PF3CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PF3CL2);
   annotation (preferredView = "info");
  end PF3CL2;

  record PF4CL "PF4CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PF4CL);
   annotation (preferredView = "info");
  end PF4CL;

  record PF5 "PF5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PF5);
   annotation (preferredView = "info");
  end PF5;

  record PH "PH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PH);
   annotation (preferredView = "info");
  end PH;

  record PH2 "PH2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PH2);
   annotation (preferredView = "info");
  end PH2;

  record PH2minus "PH2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PH2minus,
    z=-1);
   annotation (preferredView = "info");
  end PH2minus;

  record PH3 "PH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PH3);
   annotation (preferredView = "info");
  end PH3;

  record PN "PN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PN);
   annotation (preferredView = "info");
  end PN;

  record PO "PO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PO);
   annotation (preferredView = "info");
  end PO;

  record POminus "POminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.POminus,
    z=-1);
   annotation (preferredView = "info");
  end POminus;

  record POCL3 "POCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.POCL3);
   annotation (preferredView = "info");
  end POCL3;

  record POFCL2 "POFCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.POFCL2);
   annotation (preferredView = "info");
  end POFCL2;

  record POF2CL "POF2CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.POF2CL);
   annotation (preferredView = "info");
  end POF2CL;

  record POF3 "POF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.POF3);
   annotation (preferredView = "info");
  end POF3;

  record PO2 "PO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PO2);
   annotation (preferredView = "info");
  end PO2;

  record PO2minus "PO2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PO2minus,
    z=-1);
   annotation (preferredView = "info");
  end PO2minus;

  record PS "PS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PS);
   annotation (preferredView = "info");
  end PS;

  record P2 "P2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P2);
   annotation (preferredView = "info");
  end P2;

  record P2O3 "P2O3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P2O3);
   annotation (preferredView = "info");
  end P2O3;

  record P2O4 "P2O4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P2O4);
   annotation (preferredView = "info");
  end P2O4;

  record P2O5 "P2O5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P2O5);
   annotation (preferredView = "info");
  end P2O5;

  record P3 "P3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P3);
   annotation (preferredView = "info");
  end P3;

  record P3O6 "P3O6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P3O6);
   annotation (preferredView = "info");
  end P3O6;

  record P4 "P4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P4);
   annotation (preferredView = "info");
  end P4;

  record P4O6 "P4O6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O6);
   annotation (preferredView = "info");
  end P4O6;

  record P4O7 "P4O7(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O7);
   annotation (preferredView = "info");
  end P4O7;

  record P4O8 "P4O8(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O8);
   annotation (preferredView = "info");
  end P4O8;

  record P4O9 "P4O9(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O9);
   annotation (preferredView = "info");
  end P4O9;

  record P4O10 "P4O10(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O10);
   annotation (preferredView = "info");
  end P4O10;

  record Pb "Pb(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Pb);
   annotation (preferredView = "info");
  end Pb;

  record Pbplus "Pbplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Pbplus,
    z=1);
   annotation (preferredView = "info");
  end Pbplus;

  record Pbminus "Pbminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Pbminus,
    z=-1);
   annotation (preferredView = "info");
  end Pbminus;

  record PbBr "PbBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbBr);
   annotation (preferredView = "info");
  end PbBr;

  record PbBr2 "PbBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbBr2);
   annotation (preferredView = "info");
  end PbBr2;

  record PbBr3 "PbBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbBr3);
   annotation (preferredView = "info");
  end PbBr3;

  record PbBr4 "PbBr4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbBr4);
   annotation (preferredView = "info");
  end PbBr4;

  record PbCL "PbCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbCL);
   annotation (preferredView = "info");
  end PbCL;

  record PbCL2 "PbCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbCL2);
   annotation (preferredView = "info");
  end PbCL2;

  record PbCL3 "PbCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbCL3);
   annotation (preferredView = "info");
  end PbCL3;

  record PbCL4 "PbCL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbCL4);
   annotation (preferredView = "info");
  end PbCL4;

  record PbF "PbF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbF);
   annotation (preferredView = "info");
  end PbF;

  record PbF2 "PbF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbF2);
   annotation (preferredView = "info");
  end PbF2;

  record PbF3 "PbF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbF3);
   annotation (preferredView = "info");
  end PbF3;

  record PbF4 "PbF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbF4);
   annotation (preferredView = "info");
  end PbF4;

  record PbI "PbI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbI);
   annotation (preferredView = "info");
  end PbI;

  record PbI2 "PbI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbI2);
   annotation (preferredView = "info");
  end PbI2;

  record PbI3 "PbI3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbI3);
   annotation (preferredView = "info");
  end PbI3;

  record PbI4 "PbI4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbI4);
   annotation (preferredView = "info");
  end PbI4;

  record PbO "PbO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbO);
   annotation (preferredView = "info");
  end PbO;

  record PbO2 "PbO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbO2);
   annotation (preferredView = "info");
  end PbO2;

  record PbS "PbS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbS);
   annotation (preferredView = "info");
  end PbS;

  record PbS2 "PbS2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.PbS2);
   annotation (preferredView = "info");
  end PbS2;

  record Rb "Rb(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb);
   annotation (preferredView = "info");
  end Rb;

  record Rbplus "Rbplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rbplus,
    z=1);
   annotation (preferredView = "info");
  end Rbplus;

  record Rbminus "Rbminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rbminus,
    z=-1);
   annotation (preferredView = "info");
  end Rbminus;

  record RbBO2 "RbBO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbBO2);
   annotation (preferredView = "info");
  end RbBO2;

  record RbBr "RbBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbBr);
   annotation (preferredView = "info");
  end RbBr;

  record RbCL "RbCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbCL);
   annotation (preferredView = "info");
  end RbCL;

  record RbF "RbF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbF);
   annotation (preferredView = "info");
  end RbF;

  record RbH "RbH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbH);
   annotation (preferredView = "info");
  end RbH;

  record RbI "RbI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbI);
   annotation (preferredView = "info");
  end RbI;

  record RbK "RbK(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbK);
   annotation (preferredView = "info");
  end RbK;

  record RbLi "RbLi(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbLi);
   annotation (preferredView = "info");
  end RbLi;

  record RbNO2 "RbNO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbNO2);
   annotation (preferredView = "info");
  end RbNO2;

  record RbNO3 "RbNO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbNO3);
   annotation (preferredView = "info");
  end RbNO3;

  record RbNa "RbNa(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbNa);
   annotation (preferredView = "info");
  end RbNa;

  record RbO "RbO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbO);
   annotation (preferredView = "info");
  end RbO;

  record RbOH "RbOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.RbOH);
   annotation (preferredView = "info");
  end RbOH;

  record Rb2Br2 "Rb2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2Br2);
   annotation (preferredView = "info");
  end Rb2Br2;

  record Rb2CL2 "Rb2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2CL2);
   annotation (preferredView = "info");
  end Rb2CL2;

  record Rb2F2 "Rb2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2F2);
   annotation (preferredView = "info");
  end Rb2F2;

  record Rb2I2 "Rb2I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2I2);
   annotation (preferredView = "info");
  end Rb2I2;

  record Rb2O "Rb2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2O);
   annotation (preferredView = "info");
  end Rb2O;

  record Rb2O2 "Rb2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2O2);
   annotation (preferredView = "info");
  end Rb2O2;

  record Rb2O2H2 "Rb2O2H2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2O2H2);
   annotation (preferredView = "info");
  end Rb2O2H2;

  record Rb2SO4 "Rb2SO4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2SO4);
   annotation (preferredView = "info");
  end Rb2SO4;

  record Rn "Rn(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rn);
   annotation (preferredView = "info");
  end Rn;

  record Rnplus "Rnplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Rnplus,
    z=1);
   annotation (preferredView = "info");
  end Rnplus;

  record S "S(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S);
   annotation (preferredView = "info");
  end S;

  record Splus "Splus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Splus,
    z=1);
   annotation (preferredView = "info");
  end Splus;

  record Sminus "Sminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Sminus,
    z=-1);
   annotation (preferredView = "info");
  end Sminus;

  record SCL "SCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SCL);
   annotation (preferredView = "info");
  end SCL;

  record SCL2 "SCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SCL2);
   annotation (preferredView = "info");
  end SCL2;

  record SCL2plus "SCL2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SCL2plus,
    z=1);
   annotation (preferredView = "info");
  end SCL2plus;

  record SD "SD(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SD);
   annotation (preferredView = "info");
  end SD;

  record SF "SF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF);
   annotation (preferredView = "info");
  end SF;

  record SFplus "SFplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SFplus,
    z=1);
   annotation (preferredView = "info");
  end SFplus;

  record SFminus "SFminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SFminus,
    z=-1);
   annotation (preferredView = "info");
  end SFminus;

  record SF2 "SF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF2);
   annotation (preferredView = "info");
  end SF2;

  record SF2plus "SF2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF2plus,
    z=1);
   annotation (preferredView = "info");
  end SF2plus;

  record SF2minus "SF2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF2minus,
    z=-1);
   annotation (preferredView = "info");
  end SF2minus;

  record SF3 "SF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF3);
   annotation (preferredView = "info");
  end SF3;

  record SF3plus "SF3plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF3plus,
    z=1);
   annotation (preferredView = "info");
  end SF3plus;

  record SF3minus "SF3minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF3minus,
    z=-1);
   annotation (preferredView = "info");
  end SF3minus;

  record SF4 "SF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF4);
   annotation (preferredView = "info");
  end SF4;

  record SF4plus "SF4plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF4plus,
    z=1);
   annotation (preferredView = "info");
  end SF4plus;

  record SF4minus "SF4minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF4minus,
    z=-1);
   annotation (preferredView = "info");
  end SF4minus;

  record SF5 "SF5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF5);
   annotation (preferredView = "info");
  end SF5;

  record SF5plus "SF5plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF5plus,
    z=1);
   annotation (preferredView = "info");
  end SF5plus;

  record SF5minus "SF5minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF5minus,
    z=-1);
   annotation (preferredView = "info");
  end SF5minus;

  record SF6 "SF6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF6);
   annotation (preferredView = "info");
  end SF6;

  record SF6minus "SF6minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SF6minus,
    z=-1);
   annotation (preferredView = "info");
  end SF6minus;

  record SH "SH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SH);
   annotation (preferredView = "info");
  end SH;

  record SHminus "SHminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SHminus,
    z=-1);
   annotation (preferredView = "info");
  end SHminus;

  record SN "SN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SN);
   annotation (preferredView = "info");
  end SN;

  record SO "SO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SO);
   annotation (preferredView = "info");
  end SO;

  record SOminus "SOminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SOminus,
    z=-1);
   annotation (preferredView = "info");
  end SOminus;

  record SOF2 "SOF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SOF2);
   annotation (preferredView = "info");
  end SOF2;

  record SO2 "SO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2);
   annotation (preferredView = "info");
  end SO2;

  record SO2minus "SO2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2minus,
    z=-1);
   annotation (preferredView = "info");
  end SO2minus;

  record SO2CL2 "SO2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2CL2);
   annotation (preferredView = "info");
  end SO2CL2;

  record SO2FCL "SO2FCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2FCL);
   annotation (preferredView = "info");
  end SO2FCL;

  record SO2F2 "SO2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2F2);
   annotation (preferredView = "info");
  end SO2F2;

  record SO3 "SO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SO3);
   annotation (preferredView = "info");
  end SO3;

  record S2 "S2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S2);
   annotation (preferredView = "info");
  end S2;

  record S2minus "S2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S2minus,
    z=-1);
   annotation (preferredView = "info");
  end S2minus;

  record S2CL2 "S2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S2CL2);
   annotation (preferredView = "info");
  end S2CL2;

  record S2F2 "S2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S2F2);
   annotation (preferredView = "info");
  end S2F2;

  record S2O "S2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S2O);
   annotation (preferredView = "info");
  end S2O;

  record S3 "S3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S3);
   annotation (preferredView = "info");
  end S3;

  record S4 "S4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S4);
   annotation (preferredView = "info");
  end S4;

  record S5 "S5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S5);
   annotation (preferredView = "info");
  end S5;

  record S6 "S6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S6);
   annotation (preferredView = "info");
  end S6;

  record S7 "S7(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S7);
   annotation (preferredView = "info");
  end S7;

  record S8 "S8(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.S8);
   annotation (preferredView = "info");
  end S8;

  record Sc "Sc(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Sc);
   annotation (preferredView = "info");
  end Sc;

  record Scplus "Scplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Scplus,
    z=1);
   annotation (preferredView = "info");
  end Scplus;

  record Scminus "Scminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Scminus,
    z=-1);
   annotation (preferredView = "info");
  end Scminus;

  record ScO "ScO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ScO);
   annotation (preferredView = "info");
  end ScO;

  record ScOplus "ScOplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ScOplus,
    z=1);
   annotation (preferredView = "info");
  end ScOplus;

  record ScO2 "ScO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ScO2);
   annotation (preferredView = "info");
  end ScO2;

  record Sc2O "Sc2O(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Sc2O);
   annotation (preferredView = "info");
  end Sc2O;

  record Sc2O2 "Sc2O2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Sc2O2);
   annotation (preferredView = "info");
  end Sc2O2;

  record Si "Si(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Si);
   annotation (preferredView = "info");
  end Si;

  record Siplus "Siplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Siplus,
    z=1);
   annotation (preferredView = "info");
  end Siplus;

  record Siminus "Siminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Siminus,
    z=-1);
   annotation (preferredView = "info");
  end Siminus;

  record SiBr "SiBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiBr);
   annotation (preferredView = "info");
  end SiBr;

  record SiBr2 "SiBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiBr2);
   annotation (preferredView = "info");
  end SiBr2;

  record SiBr3 "SiBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiBr3);
   annotation (preferredView = "info");
  end SiBr3;

  record SiBr4 "SiBr4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiBr4);
   annotation (preferredView = "info");
  end SiBr4;

  record SiC "SiC(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiC);
   annotation (preferredView = "info");
  end SiC;

  record SiC2 "SiC2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiC2);
   annotation (preferredView = "info");
  end SiC2;

  record SiCL "SiCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiCL);
   annotation (preferredView = "info");
  end SiCL;

  record SiCL2 "SiCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiCL2);
   annotation (preferredView = "info");
  end SiCL2;

  record SiCL3 "SiCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiCL3);
   annotation (preferredView = "info");
  end SiCL3;

  record SiCL4 "SiCL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiCL4);
   annotation (preferredView = "info");
  end SiCL4;

  record SiF "SiF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiF);
   annotation (preferredView = "info");
  end SiF;

  record SiFCL "SiFCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiFCL);
   annotation (preferredView = "info");
  end SiFCL;

  record SiF2 "SiF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiF2);
   annotation (preferredView = "info");
  end SiF2;

  record SiF3 "SiF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiF3);
   annotation (preferredView = "info");
  end SiF3;

  record SiF4 "SiF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiF4);
   annotation (preferredView = "info");
  end SiF4;

  record SiH "SiH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH);
   annotation (preferredView = "info");
  end SiH;

  record SiHplus "SiHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHplus,
    z=1);
   annotation (preferredView = "info");
  end SiHplus;

  record SiHBr3 "SiHBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHBr3);
   annotation (preferredView = "info");
  end SiHBr3;

  record SiHCL "SiHCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHCL);
   annotation (preferredView = "info");
  end SiHCL;

  record SiHCL3 "SiHCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHCL3);
   annotation (preferredView = "info");
  end SiHCL3;

  record SiHF "SiHF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHF);
   annotation (preferredView = "info");
  end SiHF;

  record SiHF3 "SiHF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHF3);
   annotation (preferredView = "info");
  end SiHF3;

  record SiHI3 "SiHI3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHI3);
   annotation (preferredView = "info");
  end SiHI3;

  record SiH2 "SiH2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2);
   annotation (preferredView = "info");
  end SiH2;

  record SiH2Br2 "SiH2Br2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2Br2);
   annotation (preferredView = "info");
  end SiH2Br2;

  record SiH2CL2 "SiH2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2CL2);
   annotation (preferredView = "info");
  end SiH2CL2;

  record SiH2F2 "SiH2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2F2);
   annotation (preferredView = "info");
  end SiH2F2;

  record SiH2I2 "SiH2I2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2I2);
   annotation (preferredView = "info");
  end SiH2I2;

  record SiH3 "SiH3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3);
   annotation (preferredView = "info");
  end SiH3;

  record SiH3Br "SiH3Br(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3Br);
   annotation (preferredView = "info");
  end SiH3Br;

  record SiH3CL "SiH3CL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3CL);
   annotation (preferredView = "info");
  end SiH3CL;

  record SiH3F "SiH3F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3F);
   annotation (preferredView = "info");
  end SiH3F;

  record SiH3I "SiH3I(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3I);
   annotation (preferredView = "info");
  end SiH3I;

  record SiH4 "SiH4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH4);
   annotation (preferredView = "info");
  end SiH4;

  record SiI "SiI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiI);
   annotation (preferredView = "info");
  end SiI;

  record SiI2 "SiI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiI2);
   annotation (preferredView = "info");
  end SiI2;

  record SiN "SiN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiN);
   annotation (preferredView = "info");
  end SiN;

  record SiO "SiO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiO);
   annotation (preferredView = "info");
  end SiO;

  record SiO2 "SiO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiO2);
   annotation (preferredView = "info");
  end SiO2;

  record SiS "SiS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiS);
   annotation (preferredView = "info");
  end SiS;

  record SiS2 "SiS2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SiS2);
   annotation (preferredView = "info");
  end SiS2;

  record Si2 "Si2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Si2);
   annotation (preferredView = "info");
  end Si2;

  record Si2C "Si2C(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Si2C);
   annotation (preferredView = "info");
  end Si2C;

  record Si2F6 "Si2F6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Si2F6);
   annotation (preferredView = "info");
  end Si2F6;

  record Si2N "Si2N(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Si2N);
   annotation (preferredView = "info");
  end Si2N;

  record Si3 "Si3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Si3);
   annotation (preferredView = "info");
  end Si3;

  record Sn "Sn(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Sn);
   annotation (preferredView = "info");
  end Sn;

  record Snplus "Snplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Snplus,
    z=1);
   annotation (preferredView = "info");
  end Snplus;

  record Snminus "Snminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Snminus,
    z=-1);
   annotation (preferredView = "info");
  end Snminus;

  record SnBr "SnBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnBr);
   annotation (preferredView = "info");
  end SnBr;

  record SnBr2 "SnBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnBr2);
   annotation (preferredView = "info");
  end SnBr2;

  record SnBr3 "SnBr3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnBr3);
   annotation (preferredView = "info");
  end SnBr3;

  record SnBr4 "SnBr4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnBr4);
   annotation (preferredView = "info");
  end SnBr4;

  record SnCL "SnCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnCL);
   annotation (preferredView = "info");
  end SnCL;

  record SnCL2 "SnCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnCL2);
   annotation (preferredView = "info");
  end SnCL2;

  record SnCL3 "SnCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnCL3);
   annotation (preferredView = "info");
  end SnCL3;

  record SnCL4 "SnCL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnCL4);
   annotation (preferredView = "info");
  end SnCL4;

  record SnF "SnF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnF);
   annotation (preferredView = "info");
  end SnF;

  record SnF2 "SnF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnF2);
   annotation (preferredView = "info");
  end SnF2;

  record SnF3 "SnF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnF3);
   annotation (preferredView = "info");
  end SnF3;

  record SnF4 "SnF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnF4);
   annotation (preferredView = "info");
  end SnF4;

  record SnI "SnI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnI);
   annotation (preferredView = "info");
  end SnI;

  record SnI2 "SnI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnI2);
   annotation (preferredView = "info");
  end SnI2;

  record SnI3 "SnI3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnI3);
   annotation (preferredView = "info");
  end SnI3;

  record SnI4 "SnI4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnI4);
   annotation (preferredView = "info");
  end SnI4;

  record SnO "SnO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnO);
   annotation (preferredView = "info");
  end SnO;

  record SnO2 "SnO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnO2);
   annotation (preferredView = "info");
  end SnO2;

  record SnS "SnS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnS);
   annotation (preferredView = "info");
  end SnS;

  record SnS2 "SnS2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SnS2);
   annotation (preferredView = "info");
  end SnS2;

  record Sn2 "Sn2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Sn2);
   annotation (preferredView = "info");
  end Sn2;

  record Sr "Sr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Sr);
   annotation (preferredView = "info");
  end Sr;

  record Srplus "Srplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Srplus,
    z=1);
   annotation (preferredView = "info");
  end Srplus;

  record SrBr "SrBr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrBr);
   annotation (preferredView = "info");
  end SrBr;

  record SrBr2 "SrBr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrBr2);
   annotation (preferredView = "info");
  end SrBr2;

  record SrCL "SrCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrCL);
   annotation (preferredView = "info");
  end SrCL;

  record SrCLplus "SrCLplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrCLplus,
    z=1);
   annotation (preferredView = "info");
  end SrCLplus;

  record SrCL2 "SrCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrCL2);
   annotation (preferredView = "info");
  end SrCL2;

  record SrF "SrF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrF);
   annotation (preferredView = "info");
  end SrF;

  record SrFplus "SrFplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrFplus,
    z=1);
   annotation (preferredView = "info");
  end SrFplus;

  record SrF2 "SrF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrF2);
   annotation (preferredView = "info");
  end SrF2;

  record SrH "SrH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrH);
   annotation (preferredView = "info");
  end SrH;

  record SrI "SrI(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrI);
   annotation (preferredView = "info");
  end SrI;

  record SrI2 "SrI2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrI2);
   annotation (preferredView = "info");
  end SrI2;

  record SrO "SrO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrO);
   annotation (preferredView = "info");
  end SrO;

  record SrOH "SrOH(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrOH);
   annotation (preferredView = "info");
  end SrOH;

  record SrOHplus "SrOHplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrOHplus,
    z=1);
   annotation (preferredView = "info");
  end SrOHplus;

  record Sr_OH_2 "Sr_OH_2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Sr_OH_2);
   annotation (preferredView = "info");
  end Sr_OH_2;

  record SrS "SrS(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.SrS);
   annotation (preferredView = "info");
  end SrS;

  record Sr2 "Sr2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Sr2);
   annotation (preferredView = "info");
  end Sr2;

  record Ta "Ta(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ta);
   annotation (preferredView = "info");
  end Ta;

  record Taplus "Taplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Taplus,
    z=1);
   annotation (preferredView = "info");
  end Taplus;

  record Taminus "Taminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Taminus,
    z=-1);
   annotation (preferredView = "info");
  end Taminus;

  record TaCL5 "TaCL5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TaCL5);
   annotation (preferredView = "info");
  end TaCL5;

  record TaO "TaO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TaO);
   annotation (preferredView = "info");
  end TaO;

  record TaO2 "TaO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TaO2);
   annotation (preferredView = "info");
  end TaO2;

  record Ti "Ti(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Ti);
   annotation (preferredView = "info");
  end Ti;

  record Tiplus "Tiplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Tiplus,
    z=1);
   annotation (preferredView = "info");
  end Tiplus;

  record Timinus "Timinus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Timinus,
    z=-1);
   annotation (preferredView = "info");
  end Timinus;

  record TiCL "TiCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TiCL);
   annotation (preferredView = "info");
  end TiCL;

  record TiCL2 "TiCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TiCL2);
   annotation (preferredView = "info");
  end TiCL2;

  record TiCL3 "TiCL3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TiCL3);
   annotation (preferredView = "info");
  end TiCL3;

  record TiCL4 "TiCL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TiCL4);
   annotation (preferredView = "info");
  end TiCL4;

  record TiO "TiO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TiO);
   annotation (preferredView = "info");
  end TiO;

  record TiOplus "TiOplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TiOplus,
    z=1);
   annotation (preferredView = "info");
  end TiOplus;

  record TiOCL "TiOCL(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TiOCL);
   annotation (preferredView = "info");
  end TiOCL;

  record TiOCL2 "TiOCL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TiOCL2);
   annotation (preferredView = "info");
  end TiOCL2;

  record TiO2 "TiO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.TiO2);
   annotation (preferredView = "info");
  end TiO2;

  record U "U(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.U);
   annotation (preferredView = "info");
  end U;

  record UF "UF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF);
   annotation (preferredView = "info");
  end UF;

  record UFplus "UFplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UFplus,
    z=1);
   annotation (preferredView = "info");
  end UFplus;

  record UFminus "UFminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UFminus,
    z=-1);
   annotation (preferredView = "info");
  end UFminus;

  record UF2 "UF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF2);
   annotation (preferredView = "info");
  end UF2;

  record UF2plus "UF2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF2plus,
    z=1);
   annotation (preferredView = "info");
  end UF2plus;

  record UF2minus "UF2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF2minus,
    z=-1);
   annotation (preferredView = "info");
  end UF2minus;

  record UF3 "UF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF3);
   annotation (preferredView = "info");
  end UF3;

  record UF3plus "UF3plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF3plus,
    z=1);
   annotation (preferredView = "info");
  end UF3plus;

  record UF3minus "UF3minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF3minus,
    z=-1);
   annotation (preferredView = "info");
  end UF3minus;

  record UF4 "UF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF4);
   annotation (preferredView = "info");
  end UF4;

  record UF4plus "UF4plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF4plus,
    z=1);
   annotation (preferredView = "info");
  end UF4plus;

  record UF4minus "UF4minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF4minus,
    z=-1);
   annotation (preferredView = "info");
  end UF4minus;

  record UF5 "UF5(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF5);
   annotation (preferredView = "info");
  end UF5;

  record UF5plus "UF5plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF5plus,
    z=1);
   annotation (preferredView = "info");
  end UF5plus;

  record UF5minus "UF5minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF5minus,
    z=-1);
   annotation (preferredView = "info");
  end UF5minus;

  record UF6 "UF6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF6);
   annotation (preferredView = "info");
  end UF6;

  record UF6minus "UF6minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UF6minus,
    z=-1);
   annotation (preferredView = "info");
  end UF6minus;

  record UO "UO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UO);
   annotation (preferredView = "info");
  end UO;

  record UOplus "UOplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UOplus,
    z=1);
   annotation (preferredView = "info");
  end UOplus;

  record UOF "UOF(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UOF);
   annotation (preferredView = "info");
  end UOF;

  record UOF2 "UOF2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UOF2);
   annotation (preferredView = "info");
  end UOF2;

  record UOF3 "UOF3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UOF3);
   annotation (preferredView = "info");
  end UOF3;

  record UOF4 "UOF4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UOF4);
   annotation (preferredView = "info");
  end UOF4;

  record UO2 "UO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2);
   annotation (preferredView = "info");
  end UO2;

  record UO2plus "UO2plus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2plus,
    z=1);
   annotation (preferredView = "info");
  end UO2plus;

  record UO2minus "UO2minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2minus,
    z=-1);
   annotation (preferredView = "info");
  end UO2minus;

  record UO2F "UO2F(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2F);
   annotation (preferredView = "info");
  end UO2F;

  record UO2F2 "UO2F2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2F2);
   annotation (preferredView = "info");
  end UO2F2;

  record UO3 "UO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UO3);
   annotation (preferredView = "info");
  end UO3;

  record UO3minus "UO3minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.UO3minus,
    z=-1);
   annotation (preferredView = "info");
  end UO3minus;

  record V "V(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.V);
   annotation (preferredView = "info");
  end V;

  record Vplus "Vplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Vplus,
    z=1);
   annotation (preferredView = "info");
  end Vplus;

  record Vminus "Vminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Vminus,
    z=-1);
   annotation (preferredView = "info");
  end Vminus;

  record VCL4 "VCL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.VCL4);
   annotation (preferredView = "info");
  end VCL4;

  record VN "VN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.VN);
   annotation (preferredView = "info");
  end VN;

  record VO "VO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.VO);
   annotation (preferredView = "info");
  end VO;

  record VO2 "VO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.VO2);
   annotation (preferredView = "info");
  end VO2;

  record V4O10 "V4O10(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.V4O10);
   annotation (preferredView = "info");
  end V4O10;

  record W "W(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.W);
   annotation (preferredView = "info");
  end W;

  record Wplus "Wplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Wplus,
    z=1);
   annotation (preferredView = "info");
  end Wplus;

  record Wminus "Wminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Wminus,
    z=-1);
   annotation (preferredView = "info");
  end Wminus;

  record WCL6 "WCL6(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.WCL6);
   annotation (preferredView = "info");
  end WCL6;

  record WO "WO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.WO);
   annotation (preferredView = "info");
  end WO;

  record WOCL4 "WOCL4(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.WOCL4);
   annotation (preferredView = "info");
  end WOCL4;

  record WO2 "WO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.WO2);
   annotation (preferredView = "info");
  end WO2;

  record WO2CL2 "WO2CL2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.WO2CL2);
   annotation (preferredView = "info");
  end WO2CL2;

  record WO3 "WO3(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.WO3);
   annotation (preferredView = "info");
  end WO3;

  record WO3minus "WO3minus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.WO3minus,
    z=-1);
   annotation (preferredView = "info");
  end WO3minus;

  record Xe "Xe(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Xe);
   annotation (preferredView = "info");
  end Xe;

  record Xeplus "Xeplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Xeplus,
    z=1);
   annotation (preferredView = "info");
  end Xeplus;

  record Zn "Zn(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Zn);
   annotation (preferredView = "info");
  end Zn;

  record Znplus "Znplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Znplus,
    z=1);
   annotation (preferredView = "info");
  end Znplus;

  record Zr "Zr(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Zr);
   annotation (preferredView = "info");
  end Zr;

  record Zrplus "Zrplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Zrplus,
    z=1);
   annotation (preferredView = "info");
  end Zrplus;

  record Zrminus "Zrminus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.Zrminus,
    z=-1);
   annotation (preferredView = "info");
  end Zrminus;

  record ZrN "ZrN(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ZrN);
   annotation (preferredView = "info");
  end ZrN;

  record ZrO "ZrO(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ZrO);
   annotation (preferredView = "info");
  end ZrO;

  record ZrOplus "ZrOplus(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ZrOplus,
    z=1);
   annotation (preferredView = "info");
  end ZrOplus;

  record ZrO2 "ZrO2(g) MSL"
   extends Chemical.Interfaces.IdealGasMSL.SubstanceData(
    data = Modelica.Media.IdealGases.Common.SingleGasesData.ZrO2);
   annotation (preferredView = "info");
  end ZrO2;
  end IdealGasesMSL;
    extends Modelica.Icons.Package;

  record Silver_solid "Ag(s)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.1078682,
      z=0,
      DfH=0,
      DfG=0,
      Cp=25.4,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"});
        annotation (preferredView = "info");
  end Silver_solid;

  record Silver_aqueous "Ag+(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.1078682,
      z=1,
      DfH=105900,
      DfG=77100,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end Silver_aqueous;

  record SilverChloride_solid "AgCl(s)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.14332,
      z=0,
      DfH=-127030,
      DfG=-109720,
      Cp=50.8,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end SilverChloride_solid;

  record Calcium_aqueous "Ca++(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.0401,
      z=2,
      DfH=-542960,
      DfG=-542960 - 298.15*(33.67),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end Calcium_aqueous;

  record Chloride_aqueous "Cl-(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.03545,
      z=-1,
      DfH=-167460,
      DfG=-131170,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Chloride_aqueous;

  record CarbonMonoxide_gas "CO(g)"
   extends Chemical.Interfaces.IdealGas.SubstanceData(
      MolarWeight=0.02801,
      DfH=-110500,
      DfG=-137300,
      Cp=29.13,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.engineeringtoolbox.com/carbon-monoxide-d_975.html"});
        annotation (preferredView = "info");
  end CarbonMonoxide_gas;

  record CarbonMonoxide_aqueous "CO(aq*)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.02801,
      DfH=-276900,
      DfG=-110200,
      References={
          "Calculated from gas phase using Henry's coefficient from http://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Mask=10"});
    annotation (preferredView = "info");
  end CarbonMonoxide_aqueous;
            //  DfG = -8.314*298.15*log(0.00099/55.508)  +  -137300

  record CarbonDioxide_gas "CO2(g)"
   extends Chemical.Interfaces.IdealGas.SubstanceData(
      MolarWeight=0.044,
      DfH=-393500,
      DfG=-394400,
      Cp=37.1,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end CarbonDioxide_gas;

  record CarbonDioxide_aqueous "CO2(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      gamma=1.17385,
      MolarWeight=0.044,
      DfH=-412900,
      DfG=-386200,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end CarbonDioxide_aqueous;

  record Carbonate_aqueous "CO3--(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.06001,
      z=-2,
      DfH=-676300,
      DfG=-676300 - 298.15*(-497.065),
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Carbonate_aqueous;

  record Electrone_solid "e-(s)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=5.4857990946e-7,
      z=-1,
      DfH=0,
      DfG=0,
      Cp=0,
      References={
          "http://physics.nist.gov/cgi-bin/cuu/Value?mme, To solve standard electo-chemical cell potentials"});
        annotation (preferredView = "info");
  end Electrone_solid;

  record Iron2_aqueous "Fe++(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.05585,
      z=2,
      DfH=-87860,
      DfG=-87860 - 298.15*(-9.93),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end Iron2_aqueous;

  record Iron3_aqueous "Fe+++(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.05585,
      z=3,
      DfH=-47700,
      DfG=-47700 - 298.15*(-124.77),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end Iron3_aqueous;

  record Glucose_solid "Glu(s)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.1806,
      DfH=-1274500,
      DfG=-1274500 - 298.15*(-1220.66),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end Glucose_solid;

  record Hydrogen_gas "H2(g)"
   extends Chemical.Interfaces.IdealGas.SubstanceData(
      MolarWeight=0.00201588,
      z=0,
      DfH=0,
      DfG=0,
      Cp=28.8,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Hydrogen_gas;

  record CarbonicAcid_aqueous "H2CO3(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.062027,
      DfH=-699700,
      DfG=-699700 - 298.15*(-256.582),
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end CarbonicAcid_aqueous;

  record Water_gas "H2O(g)"
   extends Chemical.Interfaces.IdealGas.SubstanceData(
      MolarWeight=0.018015,
      DfH=-241830,
      DfG=-228590,
      Cp=33.6,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Water_gas;

  record Water_liquid_without_selfClustering "H2O(l) without self-clustering"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.018015,
      DfH=-285840,
      DfG=-237190,
      Cp=75.3,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});

      annotation (preferredView = "info", Documentation(info="<html>
<p><br><span style=\"font-family: Courier New;\">&nbsp;&nbsp;&nbsp;&nbsp;</span></p>
</html>"));
  end Water_liquid_without_selfClustering;
  //   Cv=74.539,
  // Enthalpy as in H2O(l) = with assumption that hydrogen bonds do not have significant enthaplies

  record Water_liquid "H2O(l) with self-clustering"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.018015,
      DfH=-285830,
      DfG=-227230,
      Cp=75.3,
      SelfClustering = true,
      SelfClustering_dH = -81.6348,
      SelfClustering_dS = 32.845554,
        References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
                                    //-77.95928,

  /*  SelfClustering_dH = -81.6348,
    SelfClustering_dS = 32.8344,
*/

    // S=(0 + Modelica.Constants.R*(273.15+25)*log(55.345/0.95-1))/(273.15+25),
    // SelfClustering_dS = (SelfClustering_dH + Modelica.Constants.R*(273.15+25)*log((55.345-1)/1))/(273.15+25),
    annotation (preferredView = "info", Documentation(info="<html>
<p><span style=\"font-family: Courier New;\">Even the tabulated formation Gibbs energy is DfG=-237190 there is another values because of water self-clustering. </span></p>
<p><br><span style=\"font-family: Courier New;\">New reported values for free water molecule in solution is calculated from water dissociation reaction.</span></p>
<p><span style=\"font-family: Courier New;\">&nbsp;&nbsp;&nbsp;&nbsp;</span></p>
</html>"));
  end Water_liquid;

  record Water_IceIh "H2O(s) - Ice I h"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.018015,
      DfH=-292639,
      DfG=-236590,
      Cp=37.77,
      References={"http://www1.lsbu.ac.uk/water/water_properties.html#pot"});
    annotation (preferredView = "info");
  end Water_IceIh;

  record DihydrogenPhosphate_aqueous "H2PO4-(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.095,
      z=-1,
      DfH=-1302480,
      DfG=-1302480 - 298.15*(-561.395),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end DihydrogenPhosphate_aqueous;

  record Hydronium_aqueous "H3O+(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.019022,
      z=1,
      DfH=-285840,
      DfG=-285840 - 298.15*(-163.17),
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Hydronium_aqueous;

  record PhosphoricAcid_aqueous "H3PO4(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.095,
      DfH=-1288000,
      DfG=-1288000 - 298.15*(-496.4),
      References={
          "https://en.wikipedia.org/wiki/Phosphoric_acid, https://www.researchgate.net/publication/6600409_Standard_thermodynamic_properties_of_H3PO4%28aq%29_over_a_wide_range_of_temperatures_and_pressures"});
        annotation (preferredView = "info");
  end PhosphoricAcid_aqueous;

  record Proton_aqueous "H+(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.001007,
      z=1,
      DfH=0,
      DfG=0,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Proton_aqueous;
                     // as hypothetical HA <-> H+ + A- simplification of H2O + HA <-> H3O+ + A-";

  record Bicarbonate_aqueous "HCO3-(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.06102,
      z=-1,
      DfH=-691100,
      DfG=-691100 - 298.15*(-348.82),
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Bicarbonate_aqueous;

  record Bicarbonate_blood "HCO3-(blood)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.06102,
      z=-1,
      DfH=-691100,
      DfG=-691100 - 298.15*(-348.82),
      gamma=0.79,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Bicarbonate_blood;

  record HydrogenPhosphate_aqueous "HPO4--(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.095,
      z=-2,
      DfH=-1298700,
      DfG=-1298700 - 298.15*(-686.232),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end HydrogenPhosphate_aqueous;

  record HydrogenSulfate_aqueous "HSO4-(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.097,
      z=-1,
      DfH=-885750,
      DfG=-752870,
      density=1800,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end HydrogenSulfate_aqueous;

  record Potassium_aqueous "K+(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.0391,
      z=1,
      DfH=-251200,
      DfG=-251200 - 298.15*(103.97),
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Potassium_aqueous;

  record Magnesium_aqueous "Mg++(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.0243,
      z=2,
      DfH=-461960,
      DfG=-461960 - 298.15*(-19.99),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.vias.org/genchem/standard_enthalpies_table.html"});
        annotation (preferredView = "info");
  end Magnesium_aqueous;

  record Sodium_aqueous "Na+(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.02299,
      z=1,
      DfH=-239660,
      DfG=-239660 - 298.15*(74.49),
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Sodium_aqueous;

  record Amonium_aqueous "NH4+(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.01804,
      z=1,
      DfH=-132800,
      DfG=-132800 - 298.15*(-178.77),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end Amonium_aqueous;

  record Oxygen_gas "O2(g)"
   extends Chemical.Interfaces.IdealGas.SubstanceData(
      MolarWeight=0.032,
      DfH=0,
      DfG=0,
      Cp=29.4,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Oxygen_gas;

  record Oxygen_gas_Shomate_298_6000 "O2(g) Shomate 298K–6000K"
   extends Chemical.Interfaces.IdealGasShomate.SubstanceData(
      MolarWeight=0.032,
      DfH=0,
      DfG=0,
      Cp=29.4,
      B=6.137261,
      C=-1.186521,
      D=0.09578,
      E=-0.219663,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, http://old.vscht.cz/fch/cz/pomucky/fchab/C.html"});
    annotation (preferredView = "info");
  end Oxygen_gas_Shomate_298_6000;

  record Oxygen_gas_Shomate_200_5000 "O2(g) Shomate 200K–5000K"
   extends Chemical.Interfaces.IdealGasShomate.SubstanceData(
      MolarWeight=0.032,
      DfH=0,
      DfG=0,
      Cp=29.4,
      B=-21.55543,
      C=2.456517,
      D=-0.16151,
      E=0.175056,
      X=44.837013,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, http://old.vscht.cz/fch/cz/pomucky/fchab/C.html"});
    annotation (preferredView = "info");
  end Oxygen_gas_Shomate_200_5000;
            //A=8.99044,

  record Oxygen_aqueous "O2(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.032,
      DfH=-11700,
      DfG=16320,
      References={
          "http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.pdf, https://books.google.cz/books?id=dr-VBAAAQBAJ&pg=PA156&lpg=PA156&dq=Gibbs+energy+of+formation++%22O2(aq)%22&source=bl&ots=09N5CxY7OD&sig=hbsTXQvX59vXBqHUjFVVIZQpHCA&hl=cs&sa=X&ei=sDQtVaeUMMaRsAHpzYHgAg&redir_esc=y#v=onepage&q=Gibbs%20energy%20of%20formation%20%20%22O2(aq)%22&f=false"});
        annotation (preferredView = "info");
  end Oxygen_aqueous;

  record Hydroxide_aqueous "OH-(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.017006,
      z=-1,
      DfH=-229940,
      DfG=-157300,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Hydroxide_aqueous;

  record Lead_solid "Pb(s)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.2072,
      z=0,
      DfH=0,
      DfG=0,
      Cp=26.4,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"});
        annotation (preferredView = "info");
  end Lead_solid;

  record LeadDioxide_solid "PbO2(s)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.2391988,
      z=0,
      DfH=-276600,
      DfG=-219000,
      Cp=64.6,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"});
        annotation (preferredView = "info");
  end LeadDioxide_solid;

  record LeadSulfate_solid "PbSO4(s)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.30326,
      z=0,
      DfH=-918400,
      DfG=-811200,
      Cp=103.2,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"});
        annotation (preferredView = "info");
  end LeadSulfate_solid;

  record Phosphate_aqueous "PO4---(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.095,
      z=-3,
      DfH=-1284070,
      DfG=-1284070 - 298.15*(-866.946),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end Phosphate_aqueous;

  record Sulphates_aqueous "SO4--(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.09607,
      z=-2,
      DfH=-907500,
      DfG=-907500 - 298.15*(-555.123),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
        annotation (preferredView = "info");
  end Sulphates_aqueous;

  record Ethanol_liquid "Ethanol C2H5OH(l)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.04607,
      z=0,
      DfH=-276980,
      DfG=-174180,
      Cp=112.4,
      density=789,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, https://en.wikipedia.org/wiki/Ethanol_(data_page)"});
        annotation (preferredView = "info");
  end Ethanol_liquid;
    //Some organic molecules: https://www.e-education.psu.edu/drupal6/files/be497b/pdf/Bioenergetics_AppA.pdf
  //Other source: http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf

  record Urea_aqueous "Urea(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.06006,
      z=0,
      DfH=-333189,
      DfG=-197150,
      References={"https://en.wikipedia.org/wiki/Urea"});
        annotation (preferredView = "info");
  end Urea_aqueous;

  record Globulins_aqueous "Glb(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=66.5,
      z=-4,
      DfH=0,
      DfG=0,
      References={"https://en.wikipedia.org/wiki/Human_serum_albumin"});
        annotation (preferredView = "info");
  end Globulins_aqueous;

  record Albumin_aqueous "Alb(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=66.5,
      z=-11.4,
      DfH=0,
      DfG=0,
      References={"https://en.wikipedia.org/wiki/Human_serum_albumin"});
        annotation (preferredView = "info");
  end Albumin_aqueous;

  record ADP3_aqueous "ADP3(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.427201,
      z=-3,
      DfH=0,
      DfG=0,
      References={"relative - designed only for ATP hydrolysis example"});
    annotation (preferredView = "info");
  end ADP3_aqueous;

  record ATP4_aqueous "ATP^4-(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.427201,
      z=-4,
      DfH = -1.0263e+6,
      DfG = -882161,
      References={"relative - designed only for ATP hydrolysis example"});
      // dle reakce: ATP + H2O <-> ADP + H2PO4-  (G=-30.5 kJ/mol. H=-20 kJ/mol)
    annotation (preferredView = "info");
  end ATP4_aqueous;

  record ATP3_aqueous "ATP^3-(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.427201,
      z=-4,
      DfH = -1.0263e+6,
      DfG = -919245,
      References={"relative - designed only for ATP hydrolysis example"});
      // dle reakce: ATP^4- + H+ <-> ATP^3-  (pKa=6.5, H=0 kJ/mol)
    annotation (preferredView = "info");
  end ATP3_aqueous;

  model OxygenGasOnTemperature
    Real cp1,cp2;
    Real H1,H2;
    Real S1,S2;
    Real T;
  equation
    T=200+time;
    cp1 = Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp(
      Oxygen_gas_Shomate_298_6000(), T);
    cp2 = Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp(
      Oxygen_gas_Shomate_200_5000(), T);
    H1 = Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(
      Oxygen_gas_Shomate_298_6000(), T);
    H2 = Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(
      Oxygen_gas_Shomate_200_5000(), T);
    S1 = Chemical.Interfaces.IdealGasShomate.molarEntropyPure(
      Oxygen_gas_Shomate_298_6000(), T);
    S2 = Chemical.Interfaces.IdealGasShomate.molarEntropyPure(
      Oxygen_gas_Shomate_200_5000(), T);
  end OxygenGasOnTemperature;

  record Nitrogen_gas "N2(g)"
     extends Chemical.Interfaces.IdealGas.SubstanceData(
        MolarWeight=0.0280134,
        DfH=0,
        DfG=0,
        Cp=29.1,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Type=JANAFG&Plot=on"});
      annotation (preferredView = "info");
  end Nitrogen_gas;

  record Methan_gas "CH4(g)"
   extends Chemical.Interfaces.IdealGas.SubstanceData(
      MolarWeight=0.01604246,
      z=0,
      DfH = -74848,
      DfG = -50794,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"});

    annotation (preferredView = "info");
  end Methan_gas;

  record Methan_aqueous "CH4(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.01604246,
      z=0,
      DfH = -88151,
      DfG = -34504,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=10#Solubility"});

    annotation (preferredView = "info");
  end Methan_aqueous;

  record AceticAcid_gas "CH3COOH(g)"
   extends Chemical.Interfaces.IdealGas.SubstanceData(
      MolarWeight=0.060052,
      z=0,
      DfH = -436071,
      DfG = -378978,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C64197&Mask=10#Solubility"});

    annotation (preferredView = "info");
  end AceticAcid_gas;

  record AceticAcid_aqueous "CH3COOH(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.060052,
      z=0,
      DfH = -488453,
      DfG = -399600,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"});

    annotation (preferredView = "info");
  end AceticAcid_aqueous;

  record Acetate_aqueous "CH3COO-(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.059052,
      z=-1,
      DfH = -488871,
      DfG = -372500,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"});
    annotation (preferredView = "info");
  end Acetate_aqueous;

  record Hydrogen_aqueous "H2(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.00201588,
      z=0,
      DfH=-4157,
      DfG=17740,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=10#Solubility"});
    annotation (preferredView = "info");
  end Hydrogen_aqueous;

  record Ethanol_gas "C2H5OH(g)"
   extends Chemical.Interfaces.IdealGas.SubstanceData(
      MolarWeight=0.04607,
      z=0,
      DfH = -235400,
      DfG = -168600,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"});

    annotation (preferredView = "info");
  end Ethanol_gas;

  record Ethanol_aqueous "C2H5OH(aq)"
   extends Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=0.04607,
      z=0,
      DfH=-290276,
      DfG=-181607,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Units=SI&Mask=10#Solubility"});
    annotation (preferredView = "info");
  end Ethanol_aqueous;
end Substances;
