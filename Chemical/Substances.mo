within Chemical;
package Substances "Definitions of some substances"
  package Gas
    // How to define new gas:
    // 1. Find it in Modelica.Media.IdealGases.Common.SingleGasesData -> mdata
    // 2. Find its free Gibbs energy of formation -> Gf
    //
    // or calculate it from known reaction:
    // e.g.
    // constant Chemical.Interfaces.Definition P = S + Chemical.Interfaces.Properties.processData(K,dH);
    // where K is molar-based dissociation constant and dH is consumed reaction heat (change of enthalphy)

   constant Chemical.Interfaces.Definition Unknown=Chemical.Interfaces.Definition(
        MM=1,
        DfH=0,
        DfG=0,
        Cp=1,
        phase=Chemical.Interfaces.Phase.Gas)
          "Unknown gas";

   constant Chemical.Interfaces.Definition CH4=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.CH4,
      Gf=-50720))
      "CH4(g)";

    constant Chemical.Interfaces.Definition CO=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.CO,
      Gf=-137170))
      "CO(g)";

    constant Chemical.Interfaces.Definition CO2=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.CO2,
      Gf=-394360))
      "CO2(g)";

    constant Chemical.Interfaces.Definition HCL=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.HCL,
      Gf=-95300))
      "HCL(g)";

    constant Chemical.Interfaces.Definition H2=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.H2,
      Gf=0))
      "H2(g)";

    constant Chemical.Interfaces.Definition H2O=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.H2O,
      Gf=-228570))
      "H2O(g)";

    constant Chemical.Interfaces.Definition H2O2=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.H2O2,
      Gf=-120350))
      "H2O2(g)";

    constant Chemical.Interfaces.Definition He=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.He,
      Gf=0))
      "He(g)";

    constant Chemical.Interfaces.Definition NO=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.NO,
      Gf=86550))
      "NO(g)";

    constant Chemical.Interfaces.Definition NO2=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.NO2,
      Gf=51310))
      "NO2(g)";

    constant Chemical.Interfaces.Definition N2=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData.N2,
      Gf=0))
      "N2(g)";

    constant Chemical.Interfaces.Definition O2=Chemical.Interfaces.Definition(
      data=Chemical.Interfaces.DataRecord(Modelica.Media.IdealGases.Common.SingleGasesData. O2,
      Gf=0))
      "O2(g)";

  end Gas;
    extends Modelica.Icons.Package;

  package Liquid
    // How to define new liquid:
    // 1. Find its molar mass -> MM
    // 2. Find its free enthalpy of formation -> DfH
    // 3. Find its free Gibbs energy of formation -> DfG
    //
    // or calculate it from known reaction:
    // e.g.
    // constant Chemical.Interfaces.Definition P = S + Chemical.Interfaces.Properties.processData(K,dH);
    // where K is molar-based dissociation constant and dH is consumed reaction heat (change of enthalphy)

     constant Chemical.Interfaces.Definition Unknown=Chemical.Interfaces.Definition(
        MM=1,
        DfH=0,
        DfG=0,
        Cp=1,
        phase=Chemical.Interfaces.Phase.Incompressible)
          "Unknown incompressible substance";

    // See article: http://dx.doi.org/10.35191/medsoft_2020_1_32_85_88
    constant Chemical.Interfaces.Definition H2O=Chemical.Interfaces.Definition(
          MM=0.01801528,
          DfH=-285830,
          DfG=-227230,
          Cp=75.3,
          phase=Chemical.Interfaces.Phase.Incompressible,
          SelfClustering=true,
          SelfClustering_dH=-81.6348,
          SelfClustering_dS=32.845554) "H2O(l)";

    constant Chemical.Interfaces.Definition H2OUnclustered=Chemical.Interfaces.Definition(
          MM=0.01801528,
          z=0,
          DfG=-237190,
          DfH=-285840,
          Cp=75.3,
          phase=Chemical.Interfaces.Phase.Incompressible,
          Vm=0.001*0.018015) "H2O(l) without clustering";

    constant Chemical.Interfaces.Definition Ethanol=Chemical.Interfaces.Definition(
          MM=0.04607,
          z=0,
          DfH=-276980,
          DfG=-174180,
          Cp=112.4,
          Vm=(1/789)*0.04607) "Ethanol(l)";
      //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, https://en.wikipedia.org/wiki/Ethanol_(data_page)



  end Liquid;

  package Solid
    // How to define new solid:
    // 1. Find its molar mass -> MM
    // 2. Find its free enthalpy of formation -> DfH
    // 3. Find its free Gibbs energy of formation -> DfG

     constant Chemical.Interfaces.Definition Ag=Chemical.Interfaces.Definition(
          MM=0.1078682,
          z=0,
          DfH=0,
          DfG=0,
          Cp=25.4) "Ag(s)";
       // http://www.vias.org/genchem/standard_enthalpies_table.html
       // http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf

    constant Chemical.Interfaces.Definition AgCl=Chemical.Interfaces.Definition(
          MM=0.14332,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=0,
          DfH=-127030,
          DfG=-109720,
          Cp=50.8) "AgCl(s)";

    constant Chemical.Interfaces.Definition e=Chemical.Interfaces.Definition(
          MM=5.4857990946e-7,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-1,
          DfH=0,
          DfG=0,
          Cp=0,
          Vm=1e-27) "e-(s) electrone";

    constant Chemical.Interfaces.Definition Glu=Chemical.Interfaces.Definition(
          MM=0.1806,
          phase=Chemical.Interfaces.Phase.Incompressible,
          DfH=-1274500,
          DfG=-1274500 - 298.15*(-1220.66)) "Glu(s)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition H2O_IceIh=Chemical.Interfaces.Definition(
          MM=0.018015,
          phase=Chemical.Interfaces.Phase.Incompressible,
          DfH=-292639,
          DfG=-236590,
          Cp=37.77) "H2O(s) - Ice I h";
    //  http://www1.lsbu.ac.uk/water/water_properties.html#pot;

    constant Chemical.Interfaces.Definition Pb=Chemical.Interfaces.Definition(
          MM=0.2072,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=0,
          DfH=0,
          DfG=0,
          Cp=26.4) "Pb(s)";

    constant Chemical.Interfaces.Definition PbO2=Chemical.Interfaces.Definition(
          MM=0.2391988,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=0,
          DfH=-276600,
          DfG=-219000,
          Cp=64.6) "PbO2(s)";

    constant Chemical.Interfaces.Definition PbSO4=Chemical.Interfaces.Definition(
          MM=0.30326,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=0,
          DfH=-918400,
          DfG=-811200,
          Cp=103.2) "PbSO4(s)";
  end Solid;

  package Aqueous
    // How to define new aquaeous substance:
    // 1. Find its molar mass -> MM
    // 2. Find its free enthalpy of formation -> DfH
    // 3. Find its free Gibbs energy of formation -> DfG
    //
    // or calculate it from gasseous phase and Henry's coefficient:
    // e.g.
    // constant Chemical.Interfaces.Definition O2 = Gas.O2 + Chemical.Interfaces.Properties.processData(0.0013*1/0.94,-1500*Modelica.Constants.R);
    // where 0.0013 is Henry's coefficient and 1500 K is its temperature dependence coefficient

    constant Chemical.Interfaces.Definition Agplus=Chemical.Interfaces.Definition(
          MM=0.1078682,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=1,
          DfH=105900,
          DfG=77100) "Ag+(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition Caplus2=Chemical.Interfaces.Definition(
          MM=0.0401,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=2,
          DfH=-542960,
          DfG=-542960 - 298.15*(33.67)) "Ca++(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition Clminus=Chemical.Interfaces.Definition(
          MM=0.03545,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-1,
          DfH=-167460,
          DfG=-131170) "Cl-(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition CO=Chemical.Interfaces.Definition(
          MM=0.02801,
          phase=Chemical.Interfaces.Phase.Incompressible,
          DfH=-276900,
          DfG=-110200) "CO(aq)";
  //  Calculated from gas phase using Henry's coefficient from http://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Mask=10;

    constant Chemical.Interfaces.Definition CO2=Chemical.Interfaces.Definition(
          MM=0.044,
          phase=Chemical.Interfaces.Phase.Incompressible,
          DfH=-412900,
          DfG=-386200) "CO2(aq)";
   //  http://www.vias.org/genchem/standard_enthalpies_table.html

    constant Chemical.Interfaces.Definition CO3minus2=Chemical.Interfaces.Definition(
          MM=0.06001,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-2,
          DfH=-676300,
          DfG=-676300 - 298.15*(-497.065)) "CO3--(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition Feplus2=Chemical.Interfaces.Definition(
          MM=0.05585,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=2,
          DfH=-87860,
          DfG=-87860 - 298.15*(-9.93)) "Fe++(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition Feplus3=Chemical.Interfaces.Definition(
          MM=0.05585,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=3,
          DfH=-47700,
          DfG=-47700 - 298.15*(-124.77)) "Fe+++(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition H2CO3=Chemical.Interfaces.Definition(
          MM=0.062027,
          phase=Chemical.Interfaces.Phase.Incompressible,
          DfH=-699700,
          DfG=-699700 - 298.15*(-256.582)) "H2CO3(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition H2PO4minus=Chemical.Interfaces.Definition(
          MM=0.095,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-1,
          DfH=-1302480,
          DfG=-1302480 - 298.15*(-561.395)) "H2PO4-(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition H3Oplus=Chemical.Interfaces.Definition(
          MM=0.019022,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=1,
          DfH=-285840,
          DfG=-285840 - 298.15*(-163.17)) "H3O+(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition H3PO4=Chemical.Interfaces.Definition(
          MM=0.095,
          phase=Chemical.Interfaces.Phase.Incompressible,
          DfH=-1288000,
          DfG=-1288000 - 298.15*(-496.4)) "H3PO4(aq)";
  //  https://en.wikipedia.org/wiki/Phosphoric_acid, https://www.researchgate.net/publication/6600409_Standard_thermodynamic_properties_of_H3PO4%28aq%29_over_a_wide_range_of_temperatures_and_pressures;

    constant Chemical.Interfaces.Definition Hplus=Chemical.Interfaces.Definition(
          MM=0.001007,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=1,
          DfH=0,
          DfG=0) "H+(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition HCO3minus=Chemical.Interfaces.Definition(
          MM=0.06102,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-1,
          DfH=-691100,
          DfG=-691100 - 298.15*(-348.82)) "HCO3-(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition HPO4minus2=Chemical.Interfaces.Definition(
          MM=0.095,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-2,
          DfH=-1298700,
          DfG=-1298700 - 298.15*(-686.232)) "HPO4--(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition HSO4minus=Chemical.Interfaces.Definition(
          MM=0.097,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-1,
          DfH=-885750,
          DfG=-752870,
          Vm=(1/1800)*0.097) "HSO4-(aq)";

    constant Chemical.Interfaces.Definition Kplus=Chemical.Interfaces.Definition(
          MM=0.0391,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=1,
          DfH=-251200,
          DfG=-251200 - 298.15*(103.97)) "K+(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition Mgplus2=Chemical.Interfaces.Definition(
          MM=0.0243,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=2,
          DfH=-461960,
          DfG=-461960 - 298.15*(-19.99)) "Mg++(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition Naplus=Chemical.Interfaces.Definition(
          MM=0.02299,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=1,
          DfH=-239660,
          DfG=-239660 - 298.15*(74.49)) "Na+(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition NH4plus=Chemical.Interfaces.Definition(
          MM=0.01804,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=1,
          DfH=-132800,
          DfG=-132800 - 298.15*(-178.77)) "NH4+(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition O2=Chemical.Interfaces.Definition(
          MM=0.032,
          phase=Chemical.Interfaces.Phase.Incompressible,
          DfH=-11700,
          DfG=16320) "O2(aq)";
  //  http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.pdf, https://books.google.cz/books?id=dr-VBAAAQBAJ&pg=PA156&lpg=PA156&dq=Gibbs+energy+of+formation++%22O2(aq)%22&source=bl&ots=09N5CxY7OD&sig=hbsTXQvX59vXBqHUjFVVIZQpHCA&hl=cs&sa=X&ei=sDQtVaeUMMaRsAHpzYHgAg&redir_esc=y#v=onepage&q=Gibbs%20energy%20of%20formation%20%20%22O2(aq)%22&f=false;

    constant Chemical.Interfaces.Definition OHminus=Chemical.Interfaces.Definition(
          MM=0.017006,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-1,
          DfH=-229940,
          DfG=-157300) "OH-(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition PO4minus3=Chemical.Interfaces.Definition(
          MM=0.095,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-3,
          DfH=-1284070,
          DfG=-1284070 - 298.15*(-866.946)) "PO4---(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition SO4minus2=Chemical.Interfaces.Definition(
          MM=0.09607,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-2,
          DfH=-907500,
          DfG=-907500 - 298.15*(-555.123)) "SO4--(aq)";
  //  http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf;

    constant Chemical.Interfaces.Definition Urea=Chemical.Interfaces.Definition(
          MM=0.06006,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=0,
          DfH=-333189,
          DfG=-197150) "Urea(aq)";
  //  https://en.wikipedia.org/wiki/Urea;

    constant Chemical.Interfaces.Definition Glb=Chemical.Interfaces.Definition(
          MM=66.5,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-4,
          DfH=0,
          DfG=0) "Glb(aq)";
  //  https://en.wikipedia.org/wiki/Human_serum_albumin;

    constant Chemical.Interfaces.Definition Alb=Chemical.Interfaces.Definition(
          MM=66.5,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-11.4,
          DfH=0,
          DfG=0) "Alb(aq)";
  //  https://en.wikipedia.org/wiki/Human_serum_albumin;

    constant Chemical.Interfaces.Definition ADP3=Chemical.Interfaces.Definition(
          MM=0.427201,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-3,
          DfH=0,
          DfG=0) "ADP3(aq)";
  //  relative - designed only for ATP hydrolysis example;

    constant Chemical.Interfaces.Definition ATPminus4=Chemical.Interfaces.Definition(
          MM=0.427201,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-4,
          DfH=-1.0263e+6,
          DfG=-882161) "ATP^4-(aq)";
  //  relative - designed only for ATP hydrolysis example;

    constant Chemical.Interfaces.Definition ATPminus3=Chemical.Interfaces.Definition(
          MM=0.427201,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-4,
          DfH=-1.0263e+6,
          DfG=-919245) "ATP^3-(aq)";
  //  relative - designed only for ATP hydrolysis example;

    constant Chemical.Interfaces.Definition CH4=Chemical.Interfaces.Definition(
          MM=0.01604246,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=0,
          DfH=-88151,
          DfG=-34504) "CH4(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=10#Solubility;

    constant Chemical.Interfaces.Definition CH3COOH=Chemical.Interfaces.Definition(
          MM=0.060052,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=0,
          DfH=-488453,
          DfG=-399600) "CH3COOH(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition CH3COOminus=Chemical.Interfaces.Definition(
          MM=0.059052,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=-1,
          DfH=-488871,
          DfG=-372500) "CH3COO-(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html;

    constant Chemical.Interfaces.Definition H2=Chemical.Interfaces.Definition(
          MM=0.00201588,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=0,
          DfH=-4157,
          DfG=17740) "H2(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=10#Solubility;

    constant Chemical.Interfaces.Definition C2H5OH=Chemical.Interfaces.Definition(
          MM=0.04607,
          phase=Chemical.Interfaces.Phase.Incompressible,
          z=0,
          DfH=-290276,
          DfG=-181607) "C2H5OH(aq)";
  //  http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Units=SI&Mask=10#Solubility;

  end Aqueous;



end Substances;
