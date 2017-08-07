within Chemical;
package Examples "Examples that demonstrate usage of chemical library"
extends Modelica.Icons.ExamplesPackage;

  package Substances "Definitions of substances"
      extends Modelica.Icons.Package;

    constant Chemical.Interfaces.Incompressible.SubstanceData Silver_solid(
      MolarWeight=0.1078682,
      z=0,
      DfH_25degC=0,
      DfG_25degC_1bar=0,
      Cp=25.4,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"}) "Ag(s)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Silver_aqueous(
      MolarWeight=0.1078682,
      z=1,
      DfH_25degC=105900,
      DfG_25degC_1bar=77100,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Ag+(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData SilverChloride_solid(
      MolarWeight=0.14332,
      z=0,
      DfH_25degC=-127030,
      DfG_25degC_1bar=-109720,
      Cp=50.8,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "AgCl(s)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Calcium_aqueous(
      MolarWeight=0.0401,
      z=2,
      DfH_25degC=-542960,
      DfG_25degC_1bar=-542960 - 298.15*(33.67),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Ca++(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Chloride_aqueous(
      MolarWeight=0.03545,
      z=-1,
      DfH_25degC=-167460,
      DfG_25degC_1bar=-131170,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "Cl-(aq)";

    constant Chemical.Interfaces.IdealGas.SubstanceData CarbonMonoxide_gas(
      MolarWeight=0.02801,
      DfH_25degC=-110500,
      DfG_25degC_1bar=-137300,
      Cp=29.13,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.engineeringtoolbox.com/carbon-monoxide-d_975.html"}) "CO(g)";

    constant Chemical.Interfaces.Incompressible.SubstanceData CarbonMonoxide_aqueous(
      MolarWeight=0.02801,
      DfH_25degC=-276900,
      DfG_25degC_1bar=-110200,
      References={"Calculated from gas phase using Henry's coefficient from http://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Mask=10"})
      "CO(aq*)";
      //  DfG_25degC_1bar = -8.314*298.15*log(0.00099/55.508)  +  -137300

    constant Chemical.Interfaces.IdealGas.SubstanceData CarbonDioxide_gas(
      MolarWeight=0.044,
      DfH_25degC=-393500,
      DfG_25degC_1bar=-394400,
      Cp=37.1,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "CO2(g)";

    constant Chemical.Interfaces.Incompressible.SubstanceData CarbonDioxide_aqueous(
      MolarWeight=0.044,
      DfH_25degC=-412900,
      DfG_25degC_1bar=-386200,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "CO2(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Carbonate_aqueous(
      MolarWeight=0.06001,
      z=-2,
      DfH_25degC=-676300,
      DfG_25degC_1bar=-676300 - 298.15*(-497.065),
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "CO3--(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Electrone_solid(
      MolarWeight=5.4857990946e-7,
      z=-1,
      DfH_25degC=0,
      DfG_25degC_1bar=0,
      Cp=0,
      References={
          "http://physics.nist.gov/cgi-bin/cuu/Value?mme, To solve standard electo-chemical cell potentials"}) "e-(s)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Iron2_aqueous(
      MolarWeight=0.05585,
      z=2,
      DfH_25degC=-87860,
      DfG_25degC_1bar=-87860 - 298.15*(-9.93),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Fe++(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Iron3_aqueous(
      MolarWeight=0.05585,
      z=3,
      DfH_25degC=-47700,
      DfG_25degC_1bar=-47700 - 298.15*(-124.77),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Fe+++(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Glucose_solid(
      MolarWeight=0.1806,
      DfH_25degC=-1274500,
      DfG_25degC_1bar=-1274500 - 298.15*(-1220.66),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Glu(s)";

    constant Chemical.Interfaces.IdealGas.SubstanceData Hydrogen_gas(
      MolarWeight=0.00201588,
      z=0,
      DfH_25degC=0,
      DfG_25degC_1bar=0,
      Cp=28.8,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"}) "H2(g)";

    constant Chemical.Interfaces.Incompressible.SubstanceData CarbonicAcid_aqueous(
      MolarWeight=0.062027,
      DfH_25degC=-699700,
      DfG_25degC_1bar=-699700 - 298.15*(-256.582),
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H2CO3(aq)";

    constant Chemical.Interfaces.IdealGas.SubstanceData Water_gas(
      MolarWeight=0.018015,
      DfH_25degC=-241830,
      DfG_25degC_1bar=-228590,
      Cp=33.6,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H2O(g)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Water_liquid(
      MolarWeight=0.018015,
      DfH_25degC=-285830,
      DfG_25degC_1bar=-237190,
      Cp=75.3,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H2O(l)";
   //   Cv=74.539,

    constant Chemical.Interfaces.Incompressible.SubstanceData Water_IceIh(
      MolarWeight=0.018015,
      DfH_25degC=-292639,
      DfG_25degC_1bar=-236590,
      Cp=37.77,
      References={"http://www1.lsbu.ac.uk/water/water_properties.html#pot"})
      "H2O(s) - Ice I h";

    constant Chemical.Interfaces.Incompressible.SubstanceData DihydrogenPhosphate_aqueous(
      MolarWeight=0.095,
      z=-1,
      DfH_25degC=-1302480,
      DfG_25degC_1bar=-1302480 - 298.15*(-561.395),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "H2PO4-(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Hydronium_aqueous(
      MolarWeight=0.019022,
      z=1,
      DfH_25degC=-285840,
      DfG_25degC_1bar=-285840 - 298.15*(-163.17),
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H3O+(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData PhosphoricAcid_aqueous(
      MolarWeight=0.095,
      DfH_25degC=-1288000,
      DfG_25degC_1bar=-1288000 - 298.15*(-496.4),
      References={
          "https://en.wikipedia.org/wiki/Phosphoric_acid, https://www.researchgate.net/publication/6600409_Standard_thermodynamic_properties_of_H3PO4%28aq%29_over_a_wide_range_of_temperatures_and_pressures"})
      "H3PO4(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Proton_aqueous(
      MolarWeight=0.001007,
      z=1,
      DfH_25degC=0,
      DfG_25degC_1bar=0,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H+(aq)";
               // as hypothetical HA <-> H+ + A- simplification of H2O + HA <-> H3O+ + A-";

    constant Chemical.Interfaces.Incompressible.SubstanceData Bicarbonate_aqueous(
      MolarWeight=0.06102,
      z=-1,
      DfH_25degC=-691100,
      DfG_25degC_1bar=-691100 - 298.15*(-348.82),
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "HCO3-(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Bicarbonate_blood(
      MolarWeight=0.06102,
      z=-1,
      DfH_25degC=-691100,
      DfG_25degC_1bar=-691100 - 298.15*(-348.82),
      gamma=0.79,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "HCO3-(blood)";

    constant Chemical.Interfaces.Incompressible.SubstanceData HydrogenPhosphate_aqueous(
      MolarWeight=0.095,
      z=-2,
      DfH_25degC=-1298700,
      DfG_25degC_1bar=-1298700 - 298.15*(-686.232),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "HPO4--(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData HydrogenSulfate_aqueous(
      MolarWeight=0.097,
      z=-1,
      DfH_25degC=-885750,
      DfG_25degC_1bar=-752870,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "HSO4-(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Potassium_aqueous(
      MolarWeight=0.0391,
      z=1,
      DfH_25degC=-251200,
      DfG_25degC_1bar=-251200 - 298.15*(103.97),
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "K+(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Magnesium_aqueous(
      MolarWeight=0.0243,
      z=2,
      DfH_25degC=-461960,
      DfG_25degC_1bar=-461960 - 298.15*(-19.99),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "Mg++(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Sodium_aqueous(
      MolarWeight=0.02299,
      z=1,
      DfH_25degC=-239660,
      DfG_25degC_1bar=-239660 - 298.15*(74.49),
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "Na+(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Amonium_aqueous(
      MolarWeight=0.01804,
      z=1,
      DfH_25degC=-132800,
      DfG_25degC_1bar=-132800 - 298.15*(-178.77),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "NH4+(aq)";

    constant Chemical.Interfaces.IdealGas.SubstanceData Oxygen_gas(
      MolarWeight=0.032,
      DfH_25degC=0,
      DfG_25degC_1bar=0,
      Cp=29.4,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"}) "O2(g)";

    constant Chemical.Interfaces.IdealGasShomate.SubstanceData Oxygen_gas_Shomate_298_6000(
      MolarWeight=0.032,
      DfH_25degC=0,
      DfG_25degC_1bar=0,
      cp_25degC=29.4,
      B=6.137261,
      C=-1.186521,
      D=0.09578,
      E=-0.219663,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html","http://old.vscht.cz/fch/cz/pomucky/fchab/C.html"})
      "O2(g) Shomate 298K–6000K";

    constant Chemical.Interfaces.IdealGasShomate.SubstanceData Oxygen_gas_Shomate_200_5000(
      MolarWeight=0.032,
      DfH_25degC=0,
      DfG_25degC_1bar=0,
      cp_25degC=29.4,
      B=-21.55543,
      C=2.456517,
      D=-0.16151,
      E=0.175056,
      X=44.837013,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html","http://old.vscht.cz/fch/cz/pomucky/fchab/C.html"})
      "O2(g) Shomate 200K–5000K";
      //A=8.99044,

    constant Chemical.Interfaces.Incompressible.SubstanceData Oxygen_aqueous(
      MolarWeight=0.032,
      DfH_25degC=-11700,
      DfG_25degC_1bar=16320,
      References={
          "http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.pdf, https://books.google.cz/books?id=dr-VBAAAQBAJ&pg=PA156&lpg=PA156&dq=Gibbs+energy+of+formation++%22O2(aq)%22&source=bl&ots=09N5CxY7OD&sig=hbsTXQvX59vXBqHUjFVVIZQpHCA&hl=cs&sa=X&ei=sDQtVaeUMMaRsAHpzYHgAg&redir_esc=y#v=onepage&q=Gibbs%20energy%20of%20formation%20%20%22O2(aq)%22&f=false"})
      "O2(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Hydroxide_aqueous(
      MolarWeight=0.017006,
      z=-1,
      DfH_25degC=-229940,
      DfG_25degC_1bar=-157300,
      References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "OH-(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Lead_solid(
      MolarWeight=0.2072,
      z=0,
      DfH_25degC=0,
      DfG_25degC_1bar=0,
      Cp=26.4,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"}) "Pb(s)";

    constant Chemical.Interfaces.Incompressible.SubstanceData LeadDioxide_solid(
      MolarWeight=0.2391988,
      z=0,
      DfH_25degC=-276600,
      DfG_25degC_1bar=-219000,
      Cp=64.6,
      References={
          "http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"})
      "PbO2(s)";

    constant Chemical.Interfaces.Incompressible.SubstanceData LeadSulfate_solid(
      MolarWeight=0.30326,
      z=0,
      DfH_25degC=-918400,
      DfG_25degC_1bar=-811200,
      Cp=103.2,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"})
      "PbSO4(s)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Phosphate_aqueous(
      MolarWeight=0.095,
      z=-3,
      DfH_25degC=-1284070,
      DfG_25degC_1bar=-1284070 - 298.15*(-866.946),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "PO4---(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Sulphates_aqueous(
      MolarWeight=0.09607,
      z=-2,
      DfH_25degC=-907500,
      DfG_25degC_1bar=-907500 - 298.15*(-555.123),
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "SO4--(aq)";

    constant Chemical.Interfaces.Incompressible.SubstanceData Ethanol_liquid(
      MolarWeight=0.04607,
      z=0,
      DfH_25degC=-276980,
      DfG_25degC_1bar=-174180,
      Cp=112.4,
      density=789,
      References={
          "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, https://en.wikipedia.org/wiki/Ethanol_(data_page)"})
      "Ethanol C2H5OH(l)";

      //Some organic molecules: https://www.e-education.psu.edu/drupal6/files/be497b/pdf/Bioenergetics_AppA.pdf
      //Other source: http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf
    model OxygenGasOnTemperature
      Real cp1,cp2;
      Real H1,H2;
      Real S1,S2;
      Real T;
    equation
      T=200+time;
      cp1 = Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp(Chemical.Examples.Substances.Oxygen_gas_Shomate_298_6000,T);
      cp2 = Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp(Chemical.Examples.Substances.Oxygen_gas_Shomate_200_5000,T);
      H1 = Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(Chemical.Examples.Substances.Oxygen_gas_Shomate_298_6000,T);
      H2 = Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(Chemical.Examples.Substances.Oxygen_gas_Shomate_200_5000,T);
      S1 = Chemical.Interfaces.IdealGasShomate.molarEntropyPure(Chemical.Examples.Substances.Oxygen_gas_Shomate_298_6000,T);
      S2 = Chemical.Interfaces.IdealGasShomate.molarEntropyPure(Chemical.Examples.Substances.Oxygen_gas_Shomate_200_5000,T);
    end OxygenGasOnTemperature;
  end Substances;

  model SimpleReaction
    "The simple chemical reaction A<->B with equilibrium B/A = 2"
     extends Modelica.Icons.Example;

    constant Real K = 2 "Dissociation constant of the reaction";

    constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
    constant Real R = Modelica.Constants.R "Gas constant";

    Chemical.Components.Solution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    Chemical.Components.Substance A(amountOfSubstance_start=0.9)
      annotation (Placement(transformation(extent={{-52,-8},{-32,12}})));

    Chemical.Components.Reaction reaction
      annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
    Chemical.Components.Substance B(substanceData(DfG_25degC_1bar=-R*T_25degC*
            log(K)), amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{62,-8},{42,12}})));

  equation
    connect(reaction.products[1], B.port_a) annotation (Line(
        points={{10,2},{42,2}},
        color={158,66,200},
        thickness=1));
    connect(A.solution, solution.solution) annotation (Line(
        points={{-48,-8},{-48,-92},{60,-92},{60,-98}},
        color={127,127,0}));
    connect(B.solution, solution.solution) annotation (Line(points={{58,-8},{
          58,-92},{60,-92},{60,-98}},  color={127,127,0}));
    connect(A.port_a, reaction.substrates[1]) annotation (Line(
        points={{-32,2},{-10,2}},
        color={158,66,200},
        thickness=1));
    annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Simple reaction demonstrating equilibria between substance A and substance B, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that mole fraction (A.x and B.x) are always summed to 1 for the solution.</p>
</html>"),
      experiment(StopTime=0.001));
  end SimpleReaction;

  model SimpleReaction2
    "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
     extends Modelica.Icons.Example;

    constant Real Kb(unit="kg/mol") = 2
      "Molarity based dissociation constant of the reaction with one more reactant";

    constant Real Kx(unit="1") = Kb*55.508
      "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

    constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
    constant Real R = Modelica.Constants.R "Gas constant";

    Chemical.Components.Solution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    Chemical.Components.Substance A(amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
    Chemical.Components.Reaction reaction(nS=2)
      annotation (Placement(transformation(extent={{4,-8},{24,12}})));
    Chemical.Components.Substance B(amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
    Chemical.Components.Substance C(substanceData(DfG_25degC_1bar=-R*T_25degC*
            log(Kx)), amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{68,-8},{48,12}})));

  equation
    connect(reaction.products[1], C.port_a) annotation (Line(
        points={{24,2},{48,2}},
        color={158,66,200},
        thickness=1));
    connect(A.solution, solution.solution) annotation (Line(
        points={{-30,2},{-30,-90},{60,-90},{60,-98}},
        color={127,127,0}));
    connect(C.solution, solution.solution) annotation (Line(points={{64,-8},{
          66,-8},{66,-90},{60,-90},{60,-98}},  color={127,127,0}));
    connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},
          {-30,-90},{60,-90},{60,-98}},color={127,127,0}));

    connect(B.port_a, reaction.substrates[1]) annotation (Line(
        points={{-14,-14},{-10,-14},{-10,0},{4,0}},
        color={158,66,200},
        thickness=1));
    connect(A.port_a, reaction.substrates[2]) annotation (Line(
        points={{-14,12},{-10,12},{-10,4},{4,4}},
        color={158,66,200},
        thickness=1));
    annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
      experiment(StopTime=0.001));
  end SimpleReaction2;

  model HeatingOfWater "Heating of 1 kg water"
    extends Modelica.Icons.Example;

    Chemical.Components.Solution solution(
      useElectricPort=true,
      useMechanicPorts=true,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    Chemical.Components.Substance H2O(
      redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
      substanceData=Substances.Water_liquid,
      amountOfSubstance_start=55.508)
      annotation (Placement(transformation(extent={{56,-32},{76,-12}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
      annotation (Placement(transformation(extent={{-86,-72},{-66,-52}})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-70,60},{-50,80}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
  equation
    connect(H2O.solution, solution.solution) annotation (Line(
        points={{60,-32},{60,-98}},
        color={127,127,0}));
    connect(fixedHeatFlow.port, solution.heatPort) annotation (Line(
        points={{-66,-62},{-60,-62},{-60,-102}},
        color={191,0,0}));
    connect(ground.p, solution.electricPin) annotation (Line(
        points={{-60,80},{-60,100}},
        color={0,0,255}));
  connect(fixed1.flange, solution.bottom) annotation (Line(
      points={{-18,-84},{0,-84},{0,-102}},
      color={0,127,0}));
    annotation (experiment(StopTime=1),
    Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Heating of solution by one degree, using standard HeatPort from Modelica Standard Library.</p>
<p>Observe Solution.T (or H2O.Solution.T) for temperature change.</p>
</html>"));
  end HeatingOfWater;

  model HeatingOfAlcohol "Heating of 50% ethanol"
    extends Modelica.Icons.Example;

    Chemical.Components.Solution solution(
      useElectricPort=true,
      useMechanicPorts=true,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
                              /*(mass_start=0.5 + (55.508/2)*Substances.Ethanol_liquid.MolarWeight,
    volume_start=1/(0.997*0.91251))*/
    Chemical.Components.Substance H2O(
      redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
      substanceData=Substances.Water_liquid,
      amountOfSubstance_start=55.508/2)
      annotation (Placement(transformation(extent={{-46,-8},{-26,12}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
      annotation (Placement(transformation(extent={{-86,-76},{-66,-56}})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-70,62},{-50,82}})));
    Chemical.Components.Substance Ethanol(
      redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
      substanceData=Substances.Ethanol_liquid,
      amountOfSubstance_start=55.508/2)
      annotation (Placement(transformation(extent={{18,-8},{38,12}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
  equation
    connect(H2O.solution, solution.solution) annotation (Line(
        points={{-42,-8},{-42,-34},{60,-34},{60,-98}},
        color={127,127,0}));
    connect(fixedHeatFlow.port, solution.heatPort) annotation (Line(
        points={{-66,-66},{-60,-66},{-60,-102}},
        color={191,0,0}));
    connect(ground.p, solution.electricPin) annotation (Line(
        points={{-60,82},{-60,100}},
        color={0,0,255}));
  connect(solution.solution, Ethanol.solution) annotation (Line(
      points={{60,-98},{60,-34},{22,-34},{22,-8}},
      color={127,127,0}));
  connect(solution.bottom, fixed1.flange) annotation (Line(
      points={{0,-102},{0,-84},{-18,-84}},
      color={0,127,0}));
    annotation (experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Heating of solution of water and ethanol, using standard HeatPort from Modelica Standard Library.</p>
<p>Observe Solution.T (or H2O.Solution.T or Ethanol.Solution.T) for temperature change. Note, that we can heat all substances in solution at once and the results would differ from HeatingOfWater.</p>
</html>"));
  end HeatingOfAlcohol;

  model ExothermicReaction
    "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"

     extends Modelica.Icons.Example;

    parameter Modelica.SIunits.MolarEnergy ReactionEnthalpy=-55000;

    Chemical.Components.Solution thermal_isolated_solution(
      useElectricPort=true,
      useMechanicPorts=true,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));
    Chemical.Components.Substance A(amountOfSubstance_start=0.9)
      annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
    Chemical.Components.Reaction reaction
      annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
    Chemical.Components.Substance B(amountOfSubstance_start=0.1, substanceData(
          DfH_25degC=ReactionEnthalpy))
      annotation (Placement(transformation(extent={{40,-60},{20,-40}})));

    Chemical.Components.Solution solution_at_constant_temperature(
      useElectricPort=true,
      useMechanicPorts=true,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{-100,0},{98,94}})));
    Chemical.Components.Substance A1(amountOfSubstance_start=0.9)
      annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
    Chemical.Components.Reaction reaction1
      annotation (Placement(transformation(extent={{-8,40},{12,60}})));
    Chemical.Components.Substance B1(amountOfSubstance_start=0.1, substanceData(
          DfH_25degC=ReactionEnthalpy))
      annotation (Placement(transformation(extent={{40,40},{20,60}})));

    Modelica.SIunits.HeatFlowRate q
      "Heat flow to environment to reach constant temperature";
    Modelica.SIunits.Temperature t
      "Temperature if the solution is ideally thermal isolated from environment";
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-86,56},{-66,76}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=
        298.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        origin={-78,28})));
    Chemical.Components.Substance H2O(
      redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
      substanceData=Substances.Water_liquid,
      amountOfSubstance_start=55.508)
      annotation (Placement(transformation(extent={{20,4},{40,24}})));
    Chemical.Components.Substance H2O1(
      redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
      substanceData=Substances.Water_liquid,
      amountOfSubstance_start=55.508)
      annotation (Placement(transformation(extent={{20,-94},{40,-74}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-28,4},{-8,24}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed2
      annotation (Placement(transformation(extent={{-26,-96},{-6,-76}})));
    Modelica.Electrical.Analog.Basic.Ground ground1
      annotation (Placement(transformation(extent={{-78,-32},{-58,-12}})));
  equation
    q = fixedTemperature.port.Q_flow;

    t = thermal_isolated_solution.solution.T;

    connect(A.port_a, reaction.substrates[1]) annotation (Line(
        points={{-20,-50},{-8,-50}},
        color={158,66,200},
        thickness=1));
    connect(reaction.products[1], B.port_a) annotation (Line(
        points={{12,-50},{20,-50}},
        color={158,66,200},
        thickness=1));
    connect(B.solution, thermal_isolated_solution.solution) annotation (Line(
        points={{36,-60},{36,-64},{58.4,-64},{58.4,-99.06}},
        color={127,127,0}));
    connect(A.solution, thermal_isolated_solution.solution) annotation (Line(
          points={{-36,-60},{-36,-64},{58.4,-64},{58.4,-99.06}},
                                                           color={127,127,0}));
    connect(A1.port_a, reaction1.substrates[1]) annotation (Line(
        points={{-20,50},{-8,50}},
        color={158,66,200},
        thickness=1));
    connect(reaction1.products[1], B1.port_a) annotation (Line(
        points={{12,50},{20,50}},
        color={158,66,200},
        thickness=1));
    connect(B1.solution, solution_at_constant_temperature.solution) annotation (
        Line(
        points={{36,40},{36,34},{58.4,34},{58.4,0.94}},
        color={127,127,0}));
    connect(A1.solution, solution_at_constant_temperature.solution) annotation (
        Line(points={{-36,40},{-36,34},{58.4,34},{58.4,0.94}},
                                                        color={127,127,0}));
  connect(solution_at_constant_temperature.electricPin, ground.p) annotation (
     Line(
      points={{-60.4,94},{-60,94},{-60,76},{-76,76}},
      color={0,0,255}));
  connect(fixedTemperature.port, solution_at_constant_temperature.heatPort)
    annotation (Line(
      points={{-68,28},{-60.4,28},{-60.4,-0.94}},
      color={191,0,0}));
  connect(solution_at_constant_temperature.solution, H2O.solution)
    annotation (Line(
      points={{58.4,0.94},{24,0.94},{24,4}},
      color={127,127,0}));
  connect(thermal_isolated_solution.solution, H2O1.solution) annotation (Line(
      points={{58.4,-99.06},{24,-99.06},{24,-94}},
      color={127,127,0}));
  connect(solution_at_constant_temperature.bottom, fixed1.flange) annotation (
     Line(
      points={{-1,-0.94},{0,-0.94},{0,14},{-18,14}},
      color={0,127,0}));
  connect(thermal_isolated_solution.bottom, fixed2.flange) annotation (Line(
      points={{-1,-100.94},{-1,-86},{-16,-86}},
      color={0,127,0}));
  connect(thermal_isolated_solution.electricPin, ground1.p) annotation (Line(
      points={{-60.4,-6},{-68,-6},{-68,-12}},
      color={0,0,255}));
    annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Demonstration of exotermic reaction with perfect cooling (i.e. connected fixed temperature to the HeatPort) and thermally insulated (HetPort unconnected). See solution_(...).T</p>
</html>"),
      experiment(StopTime=0.001),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics));
  end ExothermicReaction;

  model HydrogenCombustion "Hydrogen combustion in piston"
    extends Modelica.Icons.Example;

    parameter Modelica.SIunits.Volume V=0.001 "Initial volume";
   // parameter Modelica.SIunits.Pressure p=100000 "Initial pressure";
    parameter Modelica.SIunits.Temperature T=298.15 "Initial temperature";

    parameter Modelica.SIunits.Area A=0.01 "Cross area of cylinder";

    //p*V=n*R*T
   // parameter Modelica.SIunits.AmountOfSubstance n=p*V/(Modelica.Constants.R*T)
   //   "Initial amount of substances in sulution";

    Chemical.Components.Solution idealGas(
      SurfaceArea=A,
      useMechanicPorts=true,
      useThermalPort=true,
      redeclare package stateOfMatter = Interfaces.IdealGas)
      annotation (Placement(transformation(extent={{-50,-56},{50,44}})));
                     // AmbientPressure=p)
    //  volume_start=V,
    Chemical.Components.Substance H2_gas(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Substances.Hydrogen_gas,
      amountOfSubstance_start(displayUnit="mmol") = 0.026)
      annotation (Placement(transformation(extent={{-40,-32},{-20,-12}})));
    Chemical.Components.Substance O2_gas(
      substanceData=Substances.Oxygen_gas,
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      amountOfSubstance_start(displayUnit="mmol") = 0.013)
      annotation (Placement(transformation(extent={{-40,4},{-20,24}})));
    Chemical.Components.Substance H2O_gas(substanceData=Substances.Water_gas,
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas)
      annotation (Placement(transformation(extent={{44,-14},{24,6}})));
    Chemical.Components.Reaction reaction(
      nS=2,
      s={2,1},
      p={2},
      kE=4e-05)
             annotation (Placement(transformation(extent={{-10,-14},{10,6}})));
    Modelica.Mechanics.Translational.Components.Spring spring(c=1e6) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={0,58})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=2)
      annotation (Placement(transformation(extent={{-40,-86},{-20,-66}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature coolerTemperature(T=298.15)
      annotation (Placement(transformation(extent={{40,-86},{20,-66}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={0,72})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-10,-72},{10,-52}})));
  equation
  connect(reaction.products[1], H2O_gas.port_a) annotation (Line(
      points={{10,-4},{24,-4}},
      color={158,66,200},
      thickness=1));
  connect(H2_gas.port_a, reaction.substrates[1]) annotation (Line(
      points={{-20,-22},{-16,-22},{-16,-6},{-10,-6}},
      color={158,66,200},
      thickness=1));
  connect(O2_gas.port_a, reaction.substrates[2]) annotation (Line(
      points={{-20,14},{-16,14},{-16,-2},{-10,-2}},
      color={158,66,200},
      thickness=1));
  connect(H2_gas.solution, idealGas.solution) annotation (Line(
      points={{-36,-32},{30,-32},{30,-55}},
      color={127,127,0}));
  connect(O2_gas.solution, idealGas.solution) annotation (Line(
      points={{-36,4},{-44,4},{-44,-32},{30,-32},{30,-55}},
      color={127,127,0}));
  connect(H2O_gas.solution, idealGas.solution) annotation (Line(
      points={{40,-14},{40,-32},{30,-32},{30,-55}},
      color={127,127,0}));
    connect(idealGas.surfaceFlange, spring.flange_a) annotation (Line(
        points={{0,44},{0,48}},
        color={0,127,0}));
    connect(idealGas.heatPort, thermalConductor.port_a) annotation (Line(
        points={{-30,-57},{-30,-62},{-48,-62},{-48,-76},{-40,-76}},
        color={191,0,0}));
    connect(thermalConductor.port_b, coolerTemperature.port) annotation (Line(
        points={{-20,-76},{20,-76}},
        color={191,0,0}));
    connect(fixed.flange, spring.flange_b) annotation (Line(
        points={{0,72},{0,68}},
        color={0,127,0}));
  connect(idealGas.bottom, fixed1.flange) annotation (Line(
      points={{0,-57},{0,-62}},
      color={0,127,0}));
    annotation ( experiment(StopTime=1), Documentation(info="<html>
<p>The gaseous reaction of hydrogen combustion: </p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p align=\"center\"><b>2 H<sub>2</sub> + O<sub>2</sub> &lt;-&gt; 2 H<sub>2</sub>O</b></p></td>
<td><p>(1)</p></td>
</tr>
</table>
<p><br>This reaction generates a large amount of energy which can be used for mechanical or thermal purposes. </p>
<p>Building this model using the Chemical library components is easy. First, we drag and drop the library class &lsquo;Components.Solution&rsquo; into the diagram of our new model, labeled &lsquo;idealGas&rsquo; in Figure 4. In parameter dialog of this solution we check &ldquo;useThermalPorts&rdquo; and &ldquo;useMechanicsPorts&rdquo; to enable the thermal and mechanical interface. In the same dialog we need to set the area of the piston (e.g., 1 dm<sup>2</sup>), where the pressure provides the force of the green mechanical port of the uppermost side. The next parameter is the ambient external pressure surrounding the system (e.g., 1 bar). All three chemical substances of the reaction (1) can be added by dragging and dropping the library class &lsquo;Components.Substance&rsquo;. Because this model uses gases, the state of matter must be changed to some gas, such as the ideal gas prepared as &lsquo;Interfaces.IdealGas&rsquo;. The substance data must be selected to define the appropriate substances such as &lsquo;Hydrogen_gas&rsquo;, &lsquo;.Oxygen_gas&rsquo; and &lsquo;.Water_gas&rsquo; in package &lsquo;Examples.Substances&rsquo;. In addition, the initial amounts of substances can be prepared for the ideal solution of hydrogen and oxygen gases at a ratio 2:1 to attain the chemical equation above, with the expectation that at the end of the burning process, only water vapor would be presented. Therefore, the initial values of H<sub>2</sub> particles could be set to 26 mmol and of O<sub>2</sub> particles as 13 mmol. All substances must be connected with the &lsquo;idealGas&rsquo; using the blue colored solution port situated on the bottom side of each substance and solution. Then, the chemical reaction is inserted into the diagram of this model as library class &lsquo;Components.Reaction&rsquo;, and it is set to two substrates (nS=2) with stoichiometry s={2,1} and one product with stoichiometry p={2} to represent the reaction (3). The substances are then connected using violet colored substance connectors with appropriate indexes: H<sub>2</sub> to substrates[1], O<sub>2</sub> to substrates[2] and H<sub>2</sub>O to products[1]. At this point, the model is prepared to simulate the conditions of an unconnected heat port and an unconnected mechanical port. This simulation reaches the theoretical ideal of thermally isolated (zero heat flow from/to the solution) and isobaric (zero force generated on piston) conditions. </p>
<p><br><img src=\"modelica://Chemical/Resources/Images/Examples/HydrogenBurning.png\"/></p>
<p><font style=\"color: #222222; \">Mueller, M. A., Kim, T. J., Yetter, R. A., &amp; Dryer, F. L. (1999). Flow reactor studies and kinetic modeling of the H2/O2 reaction.&nbsp;<i>International Journal of Chemical Kinetics</i>,&nbsp;<i>31</i>(2), 113-125.</font></p>
<p><br>However, in the real world, there is always some thermal energy flow from the solution, and this cooling process can be connected using the thermal connector of the Modelica Standard Library 3.2.1. For example, the simple thermal conductor of thermal conductance 2W/K at a constant temperature environment of 25&deg;C is represented in the model. The mechanical power of the engine can be connected to the robust mechanical model. However, in our example we selected only a very strong mechanical spring with a spring constant of 10<sup>6</sup> N/m to stop the motion of the piston in order to generate the pressure. This standard spring component is situated above the solution in the model diagram. The results of this experiment are shown in Figure 1. </p>
</html>"));
  end HydrogenCombustion;

  model WaterVaporization "Evaporation of water"
     extends Modelica.Icons.Example;

     parameter Modelica.SIunits.Temperature T_start=273.15
      "Initial temperature";

    Chemical.Components.Solution liquid(
      temperature_start=T_start,
      useElectricPort=true,
      useMechanicPorts=true,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{-98,-98},{-6,-8}})));

    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,

    Chemical.Components.Substance H2O_liquid(substanceData=Substances.Water_liquid,
        amountOfSubstance_start=55.508) "Liquid water"
      annotation (Placement(transformation(extent={{-30,-64},{-50,-44}})));

    Chemical.Components.Solution gas(
      temperature_start=T_start,
      useMechanicPorts=true,
      useThermalPort=true,
      redeclare package stateOfMatter = Interfaces.IdealGas)
      annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                  /*volume_start(
        displayUnit="l") = 0.001, */
    Chemical.Components.GasSolubility gasSolubility(useWaterCorrection=false,
        KC=10)
      annotation (Placement(transformation(extent={{-98,24},{-78,44}})));
    Chemical.Components.Substance H2O_gaseuous(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Substances.Water_gas,
      amountOfSubstance_start(displayUnit="mmol") = 0.001)
      annotation (Placement(transformation(extent={{28,50},{8,70}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                         fixedTemperature
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        origin={84,8})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-80,-36},{-60,-16}})));
    Modelica.Blocks.Sources.Clock clock(offset=1*T_start)
      annotation (Placement(transformation(extent={{62,36},{82,56}})));
  Chemical.Components.Substance otherSubstances(substanceData=Substances.Water_liquid,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{0,28},{20,48}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(
      G=1e6) annotation (Placement(transformation(extent={{48,-8},{68,12}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-24,12},{-4,32}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed2
      annotation (Placement(transformation(extent={{-74,-92},{-54,-72}})));
  equation

    connect(H2O_liquid.solution, liquid.solution) annotation (Line(
        points={{-34,-64},{-34,-97.1},{-24.4,-97.1}},
        color={127,127,0}));
    connect(H2O_liquid.port_a, gasSolubility.liquid_port) annotation (Line(
        points={{-50,-54},{-88,-54},{-88,24}},
        color={158,66,200},
        thickness=1));
    connect(gas.solution, H2O_gaseuous.solution) annotation (Line(
        points={{27.6,6.9},{24,6.9},{24,50}},
        color={127,127,0}));
    connect(H2O_gaseuous.port_a, gasSolubility.gas_port) annotation (Line(
        points={{8,60},{-88,60},{-88,44}},
        color={158,66,200},
        thickness=1));
    connect(liquid.electricPin, ground.p) annotation (Line(
        points={{-79.6,-8},{-79.6,-8},{-70,-8},{-70,-16}},
        color={0,0,255}));
  connect(fixedTemperature.T, clock.y) annotation (Line(
      points={{96,8},{98,8},{98,46},{83,46}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(gas.solution, otherSubstances.solution) annotation (Line(
      points={{27.6,6.9},{4,6.9},{4,28}},
      color={127,127,0}));
  connect(fixedTemperature.port, thermalConductor.port_b) annotation (Line(
      points={{74,8},{72,8},{72,2},{68,2}},
      color={191,0,0}));
  connect(gas.heatPort, thermalConductor.port_a) annotation (Line(
      points={{-27.6,5.1},{-28,5.1},{-28,2},{48,2}},
      color={191,0,0}));
  connect(liquid.heatPort, thermalConductor.port_a) annotation (Line(
      points={{-79.6,-98.9},{-80,-98.9},{-80,-102},{-8,-102},{-8,2},{48,2}},
      color={191,0,0}));
  connect(liquid.bottom, fixed2.flange) annotation (Line(
      points={{-52,-98.9},{-52,-82},{-64,-82}},
      color={0,127,0}));
  connect(gas.bottom, fixed1.flange) annotation (Line(
      points={{0,5.1},{0,22},{-14,22}},
      color={0,127,0}));
    annotation (
      experiment(StopTime=100),
      Documentation(info="<html>
<p>Demonstraiton of water vaporization between two solutions - liquid and gaseous. The temperature is increased in time to illustrate, how the vaporization rate rises in higher temperatures. See liquid.T and liquid.Volume, compared to gas.T and gas.volume.</p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end WaterVaporization;

  model WaterSublimation "Sublimation of water"
     extends Modelica.Icons.Example;

     parameter Modelica.SIunits.Temperature T_start=273.15-50
      "Initial temperature";

    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,

    Chemical.Components.Solution gas(
      temperature_start=T_start,
      useMechanicPorts=true,
      useThermalPort=true,
      redeclare package stateOfMatter = Interfaces.IdealGas,
      BasePressure=600)
      annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                  /*volume_start(
        displayUnit="l") = 0.001, */
    Chemical.Components.Substance H2O_gaseuous(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Substances.Water_gas,
      amountOfSubstance_start(displayUnit="mmol") = 0.001)
      annotation (Placement(transformation(extent={{24,56},{4,76}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                         fixedTemperature
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        origin={84,8})));
    Modelica.Blocks.Sources.Clock clock(offset=1*T_start)
      annotation (Placement(transformation(extent={{62,36},{82,56}})));
  Chemical.Components.Substance otherSubstances(substanceData=Substances.Water_liquid,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-4,36},{16,56}})));
    Chemical.Components.Solution solid(
      temperature_start=T_start,
      BasePressure=600,
      useElectricPort=true,
      useMechanicPorts=true,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{8,-98},{100,-8}})));
    Chemical.Components.Substance H2O_solid(amountOfSubstance_start=55.508,
        substanceData=Substances.Water_IceIh) "Solid water"
      annotation (Placement(transformation(extent={{70,-62},{50,-42}})));
    Chemical.Components.GasSolubility gasSolubility1(useWaterCorrection=false,
        KC=10)
      annotation (Placement(transformation(extent={{-76,18},{-56,38}})));
    Modelica.Electrical.Analog.Basic.Ground ground1
      annotation (Placement(transformation(extent={{12,-38},{32,-18}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=1e6)
      annotation (Placement(transformation(extent={{48,-8},{68,12}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-18,14},{2,34}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed2
      annotation (Placement(transformation(extent={{32,-92},{52,-72}})));
  equation

    connect(gas.solution, H2O_gaseuous.solution) annotation (Line(
        points={{27.6,6.9},{20,6.9},{20,56}},
        color={127,127,0}));
  connect(fixedTemperature.T, clock.y) annotation (Line(
      points={{96,8},{98,8},{98,46},{83,46}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(gas.solution, otherSubstances.solution) annotation (Line(
      points={{27.6,6.9},{24,6.9},{24,6},{20,6},{20,36},{0,36}},
      color={127,127,0}));
    connect(solid.solution, H2O_solid.solution) annotation (Line(
        points={{81.6,-97.1},{60,-97.1},{60,-62},{66,-62}},
        color={127,127,0}));
    connect(H2O_gaseuous.port_a, gasSolubility1.gas_port) annotation (Line(
        points={{4,66},{-66,66},{-66,38}},
        color={158,66,200},
        thickness=1));
    connect(gasSolubility1.liquid_port, H2O_solid.port_a) annotation (Line(
        points={{-66,18},{-66,-52},{50,-52}},
        color={158,66,200},
        thickness=1));
    connect(solid.electricPin, ground1.p) annotation (Line(
        points={{26.4,-8},{22,-8},{22,-18}},
        color={0,0,255}));
    connect(fixedTemperature.port, thermalConductor.port_b) annotation (Line(
        points={{74,8},{72,8},{72,2},{68,2}},
        color={191,0,0}));
    connect(gas.heatPort, thermalConductor.port_a) annotation (Line(
        points={{-27.6,5.1},{-28,5.1},{-28,2},{48,2}},
        color={191,0,0}));
    connect(solid.heatPort, thermalConductor.port_a) annotation (Line(
        points={{26.4,-98.9},{-2,-98.9},{-2,2},{48,2}},
        color={191,0,0}));
  connect(solid.bottom, fixed2.flange) annotation (Line(
      points={{54,-98.9},{54,-82},{42,-82}},
      color={0,127,0}));
  connect(fixed1.flange, gas.bottom) annotation (Line(
      points={{-8,24},{0,24},{0,5.1}},
      color={0,127,0}));
    annotation (
      experiment(StopTime=50.01),
      Documentation(info="<html>
<p>Demonstraiton of water sublimation between two solutions - solid and gaseous. The temperature is increased in time to illustrate, how the sublimation rate rises in higher temperatures. See solid.T and solid.Volume, compared to gas.T and gas.volume. Note, that the liquid phase is omitted here.</p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end WaterSublimation;

  model GasSolubility "Dissolution of gases in liquids"
     extends Modelica.Icons.Example;

    Chemical.Components.Solution blood_plasma
      annotation (Placement(transformation(extent={{-100,-76},{-8,14}})));
                                        //(amountOfSolution_start=52.3)
    Chemical.Components.Solution red_cells
      annotation (Placement(transformation(extent={{8,-78},{102,14}})));
                                     //(amountOfSolution_start=39.7)

    Chemical.Components.GasSolubility CO2_dissolutionP
      annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,

    Chemical.Components.Substance CO2_unbound_plasma(substanceData=Substances.CarbonDioxide_aqueous)
      "Free dissolved CO2 in blood plasma"
      annotation (Placement(transformation(extent={{-90,-24},{-70,-4}})));
    Chemical.Components.GasSolubility O2_dissolutionP
      annotation (Placement(transformation(extent={{-34,44},{-14,64}})));

  Chemical.Sources.ExternalIdealGasSubstance O2_g_n1(
      substanceData=Substances.Oxygen_gas,
      PartialPressure=12665.626804425,
      TotalPressure=101325.0144354)
      annotation (Placement(transformation(extent={{22,76},{42,96}})));
    Chemical.Components.Substance O2_unbound_plasma(substanceData=Substances.Oxygen_aqueous)
      "Free dissolved O2 in blood plasma"
      annotation (Placement(transformation(extent={{-50,-26},{-30,-6}})));
    Chemical.Components.GasSolubility CO2_dissolutionE
      annotation (Placement(transformation(extent={{36,44},{56,64}})));

  Chemical.Sources.ExternalIdealGasSubstance CO2_g_n2(
      substanceData=Substances.CarbonDioxide_gas,
      PartialPressure=5332.8954966,
      TotalPressure=101325.0144354)
      annotation (Placement(transformation(extent={{-56,78},{-36,98}})));

    Chemical.Components.Substance CO2_unbound_erythrocyte(substanceData=
          Substances.CarbonDioxide_aqueous) "Free dissolved CO2 in red cells"
      annotation (Placement(transformation(extent={{18,-32},{38,-12}})));

    Chemical.Components.GasSolubility O2_dissolutionE_NIST(useWaterCorrection=
          true) annotation (Placement(transformation(extent={{78,44},{98,64}})));
    Chemical.Components.Substance O2_unbound_erythrocyte_NIST(substanceData=
          Substances.Oxygen_aqueous) "Free dissolved O2 in red cells"
      annotation (Placement(transformation(extent={{58,-32},{78,-12}})));
  Chemical.Components.Substance otherSubstances(substanceData=Substances.Water_liquid,
        amountOfSubstance_start=52.3)
      annotation (Placement(transformation(extent={{-42,-70},{-22,-50}})));
  Chemical.Components.Substance otherSubstances_erythrocytes(substanceData=
          Substances.Water_liquid, amountOfSubstance_start=38.7)
      annotation (Placement(transformation(extent={{64,-68},{84,-48}})));
  equation

  connect(CO2_g_n2.port_a, CO2_dissolutionP.gas_port) annotation (Line(
      points={{-36,88},{-26,88},{-26,72},{-68,72},{-68,64}},
      color={158,66,200},
      thickness=1));
  connect(CO2_g_n2.port_a, CO2_dissolutionE.gas_port) annotation (Line(
      points={{-36,88},{-26,88},{-26,72},{46,72},{46,64}},
      color={158,66,200},
      thickness=1));
  connect(CO2_dissolutionP.liquid_port, CO2_unbound_plasma.port_a)
    annotation (Line(
      points={{-68,44},{-68,-14},{-70,-14}},
      color={158,66,200},
      thickness=1));
  connect(CO2_dissolutionE.liquid_port, CO2_unbound_erythrocyte.port_a)
    annotation (Line(
      points={{46,44},{46,-22},{38,-22}},
      color={158,66,200},
      thickness=1));
  connect(O2_g_n1.port_a, O2_dissolutionP.gas_port) annotation (Line(
      points={{42,86},{66,86},{66,70},{-24,70},{-24,64}},
      color={158,66,200},
      thickness=1));
  connect(O2_dissolutionP.liquid_port, O2_unbound_plasma.port_a) annotation (
      Line(
      points={{-24,44},{-24,-16},{-30,-16}},
      color={158,66,200},
      thickness=1));
  connect(CO2_unbound_plasma.solution, blood_plasma.solution) annotation (
      Line(
      points={{-86,-24},{-86,-75.1},{-26.4,-75.1}},
      color={127,127,0}));
  connect(O2_unbound_plasma.solution, blood_plasma.solution) annotation (Line(
      points={{-46,-26},{-46,-75.1},{-26.4,-75.1}},
      color={127,127,0}));
  connect(CO2_unbound_erythrocyte.solution, red_cells.solution) annotation (
      Line(
      points={{22,-32},{22,-77.08},{83.2,-77.08}},
      color={127,127,0}));
  connect(O2_g_n1.port_a, O2_dissolutionE_NIST.gas_port) annotation (Line(
      points={{42,86},{66,86},{66,70},{88,70},{88,64}},
      color={158,66,200},
      thickness=1));
  connect(O2_dissolutionE_NIST.liquid_port, O2_unbound_erythrocyte_NIST.port_a)
    annotation (Line(
      points={{88,44},{88,-22},{78,-22}},
      color={158,66,200},
      thickness=1));
  connect(O2_unbound_erythrocyte_NIST.solution, red_cells.solution)
    annotation (Line(
      points={{62,-32},{62,-77.08},{83.2,-77.08}},
      color={127,127,0}));
  connect(blood_plasma.solution, otherSubstances.solution) annotation (Line(
      points={{-26.4,-75.1},{-38,-75.1},{-38,-70}},
      color={127,127,0}));
  connect(red_cells.solution, otherSubstances_erythrocytes.solution)
    annotation (Line(
      points={{83.2,-77.08},{68,-77.08},{68,-68}},
      color={127,127,0}));
    annotation (
      experiment(StopTime=1e-005),
      Documentation(info="<html>
<p>Demonstration of different blood gases solubility in erythrocytes and in plasma. The difference is governed by various amount of other substances in the solution. </p>
<p>Please note, that the total content of CO2 and O2 in blood plasma and erythrocytes must be determined by including bicarbonate and hemoglobin connected amounts. </p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end GasSolubility;

  model EnzymeKinetics "Basic enzyme kinetics"
    extends Modelica.Icons.Example;

    Chemical.Components.Solution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
    Chemical.Components.Substance P(substanceData(DfG_25degC_1bar=-Modelica.Constants.R
            *298.15*50))
      annotation (Placement(transformation(extent={{92,-12},{72,8}})));
    Chemical.Components.Substance S(amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));

       parameter Modelica.SIunits.AmountOfSubstance tE=1e-6
      "Total amount of enzyme";
       parameter Real k_cat(unit="1/s", displayUnit="1/min")= 1
      "Forward rate of second reaction";
       constant Modelica.SIunits.Concentration Km=0.1
      "Michaelis constant = substrate concentration at rate of half Vmax";

      parameter Modelica.SIunits.MolarFlowRate Vmax=1e-5*k_cat
      "Maximal molar flow";
      parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution= 55.508
      "Amount of solution used in kinetics";

        Chemical.Components.Substance ES(substanceData(DfG_25degC_1bar=-
            Modelica.Constants.R*298.15*log(2/Km)), amountOfSubstance_start=tE/
          2) annotation (Placement(transformation(extent={{-12,-10},{8,10}})));
        Chemical.Components.Substance E(amountOfSubstance_start=tE/2)
      annotation (Placement(transformation(extent={{-10,38},{10,58}})));
    Chemical.Components.Reaction chemicalReaction(nS=2, KC=Vmax/(2*Modelica.Constants.R
          *298.15*log(2)))
      annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));

    Chemical.Components.Reaction chemicalReaction1(nP=2, KC=Vmax/(2*Modelica.Constants.R
          *298.15*(50 - log(2))))
      annotation (Placement(transformation(extent={{24,-10},{44,10}})));

  Chemical.Components.Substance otherSubstances(substanceData=Substances.Water_liquid,
        amountOfSubstance_start=52.3)
      annotation (Placement(transformation(extent={{42,-76},{62,-56}})));
  equation
       //Michaelis-Menton: v=((E.q_out.conc + ES.q_out.conc)*k_cat)*S.concentration/(Km+S.concentration);

    connect(S.port_a, chemicalReaction.substrates[1]) annotation (Line(
        points={{-72,-4},{-56,-4},{-56,-2},{-42,-2}},
        color={158,66,200},
        thickness=1));
    connect(chemicalReaction.products[1], ES.port_a) annotation (Line(
        points={{-22,0},{8,0}},
        color={158,66,200},
        thickness=1));
    connect(ES.port_a, chemicalReaction1.substrates[1]) annotation (Line(
        points={{8,0},{24,0}},
        color={158,66,200},
        thickness=1));
    connect(E.port_a, chemicalReaction.substrates[2]) annotation (Line(
        points={{10,48},{-52,48},{-52,2},{-42,2}},
        color={158,66,200},
        thickness=1));
    connect(E.port_a, chemicalReaction1.products[2]) annotation (Line(
        points={{10,48},{54,48},{54,2},{44,2}},
        color={158,66,200},
        thickness=1));
    connect(chemicalReaction1.products[1], P.port_a) annotation (Line(
        points={{44,-2},{58,-2},{58,-2},{72,-2}},
        color={158,66,200},
        thickness=1));
    connect(E.solution, solution.solution) annotation (Line(
        points={{-6,38},{-8,38},{-8,-98},{60,-98}},
        color={127,127,0}));
    connect(ES.solution, solution.solution)
      annotation (Line(points={{-8,-10},{-8,-98},{60,-98}},         color={127,127,0}));

    connect(S.solution, solution.solution) annotation (Line(
        points={{-88,-14},{-88,-56},{-8,-56},{-8,-98},{60,-98}},
        color={127,127,0}));
    connect(P.solution, solution.solution) annotation (Line(
        points={{88,-12},{88,-98},{60,-98}},
        color={127,127,0}));
  connect(solution.solution, otherSubstances.solution) annotation (Line(
      points={{60,-98},{46,-98},{46,-76}},
      color={127,127,0}));
        annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Be carefull, the assumption for Michaelis-Menton are very strong: </p>
<p>The substrate must be in sufficiently high concentration and the product must be in very low concentration to reach almost all enzyme in enzyme-substrate complex all time. ([S] &gt;&gt; Km) &amp;&amp; ([P] &lt;&lt; K2)</p>
<p><br>To recalculate the enzyme kinetics from Michaelis-Menton parameters Km, tE a k_cat is selected the same half-rate of the reaction defined as:</p>
<p>E = ES = tE/2 .. the amount of free enzyme is the same as the amount of enzyme-substrate complexes</p>
<p>S = Km .. the amount of substrate is Km</p>
<p>r = Vmax/2 = tE*k_cat / 2 .. the rate of reaction is the half of maximal rate</p>
<p><br>Conversions of molar concentration to mole fraction (MM is molar mass of the solvent in solution -&gt; 55.508 kg/mol for water):</p>
<p>x(Km) = Km/MM</p>
<p>x(tE) = tE/MM</p>
<p>xS = S/MM = Km/MM</p>
<p><br>The new kinetics of the system defined as:</p>
<p>uS&deg; = DfG(S) = 0</p>
<p>uE&deg; = DfG(E) = 0</p>
<p>uES&deg; = <b>DfG(ES) = DfG(S) + DfG(E) - R*T*ln(2/x(Km))</b></p>
<p>from dissociation coeficient of the frist reaction 2/x(Km) = xSE/(xS*xE) = exp((uE&deg; + uS&deg; - uES&deg;)/(RT))</p>
<p>uP&deg; = DfG(P) </p>
<p><br>r = Vmax/2</p>
<p>r = -kC1 * (uES&deg; - uE&deg; - uS&deg; + R*T*ln(xES/(xE*xS) ) = -kC1 * (-R*T*ln(2/x(Km)) + R*T*ln(xS) ) = kC1 * R * T * ln(2)</p>
<p>because xES=xE this time</p>
<p>r = -kC2 * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) = -kC2 * (DfG(P) - uES&deg; + R*T*ln(xP) ) = kC2 * (-DfG(P) - R * T * ln(2))</p>
<h4>kC1 = (Vmax/2) / (R * T * ln(2))</h4>
<h4>kC2 = (Vmax/2) / ( -DfG(P) - R * T * ln(2) ) </h4>
<p><br>For example in case of C=AmountOfSolution/(Tau*ActivationPotential) we can rewrite C to ActivationPotential (Be carefull: this energy is not the same as in <a href=\"http://en.wikipedia.org/wiki/Arrhenius_equation\">Arrhenius equation</a> or in Transition State Theory):</p>
<p>ActivationPotential1 = AmountOfSolution/(Tau*(Vmax/2)) * R * T * ln(2) </p>
<p>ActivationPotential2 = AmountOfSolution/(Tau*(Vmax/2)) * ( -DfG(P) - R * T * ln(2) ) </p>
<p><br>where</p>
<p>AmountOfSolution = MM = 55.508 (for water)</p>
<p>Tau = 1 s (just to be physical unit correct)</p>
<p>DfG(P) = -R*T*50 is Gibbs energy of formation of product (setting negative enough makes second reaction almost irreversible)</p>
<h4>The maximum of the new enzyme kinetics</h4>
<p>The enzymatic rate must have a maximum near of Vmax. </p>
<p>The new maximum is a litle higher: Vmax * (1 + 1/( -uP&deg;/(R*T*ln(2)) - 1) ), for example if -uP&deg;/RT = 50, the new maximum is around 1.014*Vmax, where Vmax is the maximum of Michaelis Menten.</p>
<p>The proof:</p>
<p>We want to sutisfied the following inequality:</p>
<p>-kC2 * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) ?=&lt;? Vmax * (1 + 1/( -uP&deg;/(R*T*ln(2)) - 1) )</p>
<p><br>(Vmax/2) * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) / ( - uP&deg; - R * T * ln(2) ) ?=&lt;? Vmax*(1 + R*T*ln(2) / ( -uP&deg; - R*T*ln(2)) )</p>
<p>(uP&deg; +<b> </b>R*T*ln(2/x(Km)) + R*T*ln(xP*xE/xES) ) ?=&lt;? 2*( - uP&deg; - R * T * ln(2) ) + 2*R*T*ln(2)</p>
<p>R*T*ln(xP*xE/xES) ?=&lt;? - uP&deg; - R*T*ln(2/x(Km)) </p>
<p>xP*xE/xES ?=&lt;? exp((- uP&deg; - R*T*ln(2/x(Km))/(R*T))</p>
<p>The equality is the equation of the equilibrium: xP*xE/xES = exp((- uP&deg; - uE&deg; + uES&deg; )/(R*T)) = exp((- uP&deg; - R*T*ln(2/x(Km))/(R*T))</p>
<p>If the equilibrium of the reaction is reached only by forward rate then xP*xE/xES must be less than the dissociation constant.</p>
<h4>The increasing of the amount of the enzyme</h4>
<p>In the situation of doubled amount of enzyme should double also the maximal speed of the reaction, shouldn&apos;t?</p>
<p>The assumptions of</p>
</html>"),
      experiment(StopTime=199000));
  end EnzymeKinetics;

  model ElectrochemicalCell
    "The electrochemical cell: Pt(s) | H2(g) | H+(aq), Cl-(aq) | AgCl(s) | Ag(s)"
   extends Modelica.Icons.Example;

    Chemical.Components.Solution cathode(ElectricGround=false)
      annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
    Chemical.Components.Solution anode(ElectricGround=false)
      annotation (Placement(transformation(extent={{62,-50},{96,50}})));

    Chemical.Components.Solution solution1(ElectricGround=false)
      annotation (Placement(transformation(extent={{-30,-60},{38,6}})));

    Chemical.Components.Substance Ag(substanceData=Substances.Silver_solid,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-72,-30},{-52,-10}})));
    Chemical.Components.Substance Cl(substanceData=Substances.Chloride_aqueous,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-2,-36},{-22,-16}})));
    Chemical.Components.Substance AgCl(substanceData=Substances.SilverChloride_solid)
      annotation (Placement(transformation(extent={{-76,4},{-56,24}})));
  Chemical.Sources.ExternalIdealGasSubstance H2(
      substanceData=Substances.Hydrogen_gas,
      PartialPressure=100000,
      TotalPressure=100000)
      annotation (Placement(transformation(extent={{24,32},{44,52}})));
    Chemical.Components.Substance H(substanceData=Substances.Proton_aqueous,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{6,-36},{26,-16}})));
    Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
      annotation (Placement(transformation(extent={{-6,64},{14,84}})));
    Chemical.Components.Reaction electrodeReaction(nP=2, p={2,2}) annotation (
        Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={52,6})));
    Chemical.Components.Reaction electrodeReaction1(nS=2, nP=2) annotation (
        Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=90,
          origin={-40,0})));

  Chemical.Components.ElectronTransfer electrone
      annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                               //(substanceData=Chemical.Examples.Substances.Electrone_solid)
  Chemical.Components.ElectronTransfer electrone1
      annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                //(substanceData=Chemical.Examples.Substances.Electrone_solid)

  Modelica.Electrical.Analog.Basic.Ground ground
    annotation (Placement(transformation(extent={{84,-84},{104,-64}})));
  equation
    connect(Ag.port_a, electrodeReaction1.substrates[1]) annotation (Line(
        points={{-52,-20},{-42,-20},{-42,-10},{-42,-10}},
        color={158,66,200},
        thickness=1));
    connect(Cl.port_a, electrodeReaction1.substrates[2]) annotation (Line(
        points={{-22,-26},{-38,-26},{-38,-10}},
        color={158,66,200},
        thickness=1));
    connect(AgCl.port_a, electrodeReaction1.products[1]) annotation (Line(
        points={{-56,14},{-42,14},{-42,10}},
        color={158,66,200},
        thickness=1));
    connect(H2.port_a, electrodeReaction.substrates[1]) annotation (Line(
        points={{44,42},{52,42},{52,16}},
        color={158,66,200},
        thickness=1));
    connect(H.port_a, electrodeReaction.products[1]) annotation (Line(
        points={{26,-26},{54,-26},{54,-4}},
        color={158,66,200},
        thickness=1));
    connect(electrodeReaction.products[2], electrone1.port_a) annotation (Line(
        points={{50,-4},{50,-16},{68,-16}},
        color={158,66,200},
        thickness=1));
    connect(electrodeReaction1.products[2], electrone.port_a) annotation (Line(
        points={{-38,10},{-38,42},{-58,42}},
        color={158,66,200},
        thickness=1));
    connect(Cl.solution, solution1.solution) annotation (Line(
        points={{-6,-36},{-6,-40},{24.4,-40},{24.4,-59.34}},
        color={127,127,0}));
    connect(H.solution, solution1.solution) annotation (Line(points={{10,-36},
          {10,-40},{24.4,-40},{24.4,-59.34}},
                                       color={127,127,0}));
  connect(electrone.solution, cathode.solution) annotation (Line(
      points={{-74,32},{-74,-34},{-68,-34},{-68,-42.84},{-54.4,-42.84}},
      color={127,127,0}));
  connect(electrone1.solution, anode.solution) annotation (Line(
      points={{84,-26},{84,-49},{89.2,-49}},
      color={127,127,0}));
  connect(AgCl.solution, cathode.solution) annotation (Line(
      points={{-72,4},{-74,4},{-74,-34},{-68,-34},{-68,-42.84},{-54.4,-42.84}},
      color={127,127,0}));
  connect(Ag.solution, cathode.solution) annotation (Line(
      points={{-68,-30},{-68,-42.84},{-54.4,-42.84}},
      color={158,66,200}));
    connect(voltageSensor.p, electrone.pin) annotation (Line(
        points={{-6,74},{-96,74},{-96,42},{-78,42}},
        color={0,0,255}));
    connect(voltageSensor.n, electrone1.pin) annotation (Line(
        points={{14,74},{92,74},{92,-16},{88,-16}},
        color={0,0,255}));
    connect(electrone1.pin, ground.p) annotation (Line(
        points={{88,-16},{92,-16},{92,-64},{94,-64}},
        color={0,0,255}));
    annotation (
    experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ElectrochemicalCell;

  model RedCellMembrane
   // import Chemical;
    extends Modelica.Icons.Example;

    parameter Real KC=1;//e-6 "Slow down factor";

    Chemical.Components.Solution blood_erythrocytes(ElectricGround=false)
      annotation (Placement(transformation(extent={{-180,-100},{180,-10}})));
    Chemical.Components.Solution blood_plasma
      annotation (Placement(transformation(extent={{-180,12},{180,100}})));

    Chemical.Components.Substance HCO3(substanceData=Substances.Bicarbonate_blood,
        amountOfSubstance_start(displayUnit="mmol") = 0.024) annotation (
        Placement(transformation(extent={{-10,-10},{10,10}}, origin={-18,30})));

    Chemical.Components.Substance H2O(substanceData=Substances.Water_liquid,
        amountOfSubstance_start=51.8*0.994648)
      annotation (Placement(transformation(extent={{-146,20},{-166,40}})));
    Chemical.Components.Substance HCO3_E(substanceData=Substances.Bicarbonate_blood,
        amountOfSubstance_start(displayUnit="mmol") = 0.0116)
      annotation (Placement(transformation(extent={{-28,-38},{-8,-18}})));
    Chemical.Components.Substance H2O_E(substanceData=Substances.Water_liquid,
        amountOfSubstance_start=38.7*0.994648)
      annotation (Placement(transformation(extent={{-144,-38},{-164,-18}})));
    Chemical.Components.Substance Cl_E(substanceData=Substances.Chloride_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.0499)
      annotation (Placement(transformation(extent={{-4,-38},{16,-18}})));
    Chemical.Components.Substance Cl(substanceData=Substances.Chloride_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.123)
      annotation (Placement(transformation(extent={{-4,20},{16,40}})));

  //  Real pH_e; //,pH_p;

    Chemical.Components.Membrane Aquapirin(KC=KC) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-168,0})));
    Chemical.Components.Membrane Band3(KC=KC) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-6,0})));
    Chemical.Components.Membrane Band3_(useKineticsInput=false, KC=KC)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={18,0})));
    Chemical.Components.Substance permeableUncharged(amountOfSubstance_start(
          displayUnit="mmol") = 0.0118)
      annotation (Placement(transformation(extent={{166,20},{146,40}})));
    Chemical.Components.Substance permeableUncharged_E(amountOfSubstance_start(
          displayUnit="mmol") = 0.00903, substanceData(MolarWeight=0.1))
      annotation (Placement(transformation(extent={{164,-38},{144,-18}})));
    Chemical.Components.Substance chargedImpermeable_E(amountOfSubstance_start(
          displayUnit="mmol") = 0.0165, substanceData(MolarWeight=1))
      annotation (Placement(transformation(extent={{144,-62},{164,-42}})));
    Chemical.Components.Membrane leak(useKineticsInput=false, KC=KC)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={140,0})));
    Chemical.Components.Substance Lac_E(substanceData=Substances.Chloride_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.00062)
      annotation (Placement(transformation(extent={{56,-38},{76,-18}})));
    Chemical.Components.Substance Lac(substanceData=Substances.Chloride_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.00131)
      annotation (Placement(transformation(extent={{56,20},{76,40}})));
    Chemical.Components.Membrane MCT_(useKineticsInput=false, KC=KC)
      "Monocarboxylate transporters" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={78,0})));
    Chemical.Components.Substance H_E(substanceData=Substances.Proton_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 38.7*10^(-7.2)) "H+"
      annotation (Placement(transformation(extent={{30,-38},{50,-18}})));
    Chemical.Components.Substance H(substanceData=Substances.Proton_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 51.8*10^(-7.4))
      "H+ in plasma"
      annotation (Placement(transformation(extent={{30,20},{50,40}})));
    Chemical.Components.Membrane MCT(useKineticsInput=false, KC=KC)
      "Monocarboxylate transporters" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={52,0})));
    Chemical.Components.Substance CO2(substanceData=Substances.CarbonDioxide_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.00167)
      "free dissolved unbound CO2"
      annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
    Chemical.Components.Substance CO2_E(substanceData=Substances.CarbonDioxide_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.00125)
      "free dissolved unbound CO2"
      annotation (Placement(transformation(extent={{-58,-38},{-38,-18}})));
    Chemical.Components.Membrane freeCO2(KC=KC) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-36,0})));
    Chemical.Components.Substance O2(substanceData=Substances.Oxygen_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.000167)
      "free dissolved undound oxygen"
      annotation (Placement(transformation(extent={{96,20},{116,40}})));
    Chemical.Components.Membrane freeO2(KC=KC) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={118,0})));
    Chemical.Components.Substance O2_E(amountOfSubstance_start(displayUnit=
            "mmol") = 0.000125, substanceData=Substances.Oxygen_aqueous)
      "free dissolved undound O2"
      annotation (Placement(transformation(extent={{96,-38},{116,-18}})));
    Chemical.Components.Substance
                         K(substanceData=Substances.Potassium_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.004)
      annotation (Placement(transformation(extent={{-92,20},{-112,40}})));
    Chemical.Components.Substance
                         Na(                               substanceData=
          Substances.Sodium_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.138)
      annotation (Placement(transformation(extent={{-118,20},{-138,40}})));
    Chemical.Components.Substance
                         Na_E(substanceData=Substances.Sodium_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.007)
      annotation (Placement(transformation(extent={{-118,-38},{-138,-18}})));
    Chemical.Components.Substance
                         K_E(substanceData=Substances.Potassium_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.096)
      annotation (Placement(transformation(extent={{-112,-38},{-92,-18}})));
    Chemical.Components.Substance H2PO4_E(substanceData=Substances.DihydrogenPhosphate_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.000175)
      annotation (Placement(transformation(extent={{-84,-38},{-64,-18}})));
    Chemical.Components.Substance ADP_E(substanceData(z=-3),
        amountOfSubstance_start(displayUnit="mmol") = 9.6e-05)
      annotation (Placement(transformation(extent={{-114,-62},{-94,-42}})));
    Chemical.Components.Substance ATP_E(substanceData(
        z=-4,
        DfH_25degC=16700,
        DfG_25degC_1bar=30500,
        References={"http://www.wiley.com/college/pratt/0471393878/student/review/thermodynamics/7_relationship.html"}),
        amountOfSubstance_start(displayUnit="mmol") = 0.00128)
      annotation (Placement(transformation(extent={{-118,-62},{-138,-42}})));
    Chemical.Components.Substance HPO4_E(substanceData=Substances.HydrogenPhosphate_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.000495)
      annotation (Placement(transformation(extent={{-84,-62},{-64,-42}})));
    Chemical.Components.Substance H2PO4(substanceData=Substances.DihydrogenPhosphate_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.000365)
      annotation (Placement(transformation(extent={{-86,20},{-66,40}})));
    Chemical.Components.Substance HPO4(substanceData=Substances.HydrogenPhosphate_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.001635)
      annotation (Placement(transformation(extent={{-86,42},{-66,62}})));
    Chemical.Components.Substance albumin(substanceData(
        MolarWeight=66.463,
        z=-17,
        density=1080), amountOfSubstance_start(displayUnit="mmol") = 0.0007)
      annotation (Placement(transformation(extent={{116,76},{96,96}})));
    Chemical.Components.Substance globulins(substanceData(
        MolarWeight=34,
        z=-2.43,
        density=1080), amountOfSubstance_start(displayUnit="mmol") = 0.00082)
      annotation (Placement(transformation(extent={{150,76},{130,96}})));
    Chemical.Components.Substance Ca(substanceData=Substances.Calcium_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.00175) "Ca2+"
      annotation (Placement(transformation(extent={{-112,42},{-92,62}})));
    Chemical.Components.Substance Mg(substanceData=Substances.Magnesium_aqueous,
        amountOfSubstance_start(displayUnit="mmol") = 0.00108) "Mg2+"
      annotation (Placement(transformation(extent={{-112,-84},{-92,-64}})));
    Chemical.Components.Substance hemoglobin(substanceData(
        MolarWeight=64,
        z=-4.4,
        density=1500), amountOfSubstance_start(displayUnit="mmol") = 0.00513)
      annotation (Placement(transformation(extent={{94,-94},{74,-74}})));
    Chemical.Components.Substance DPG(amountOfSubstance_start(displayUnit=
            "mmol") = 0.0051, substanceData(
        MolarWeight=0.266,
        z=-2.2,
        density=1000))
      annotation (Placement(transformation(extent={{128,-94},{108,-74}})));
    Chemical.Components.Substance GSH(substanceData(
        MolarWeight=0.2,
        z=-1,
        density=1000), amountOfSubstance_start(displayUnit="mmol") = 0.00223)
      annotation (Placement(transformation(extent={{164,-94},{144,-74}})));
  equation
  //  pH_p = -log10(H.a);
   // pH_e = -log10(H_E.a);
  connect(H2O.solution, blood_plasma.solution)
    annotation (Line(points={{-150,20},{108,20},{108,12.88}},
                                                      color={127,127,0}));
  connect(Cl.solution, blood_plasma.solution) annotation (Line(
      points={{0,20},{0,12.88},{108,12.88}},
      color={127,127,0}));
    connect(H2O_E.solution, blood_erythrocytes.solution) annotation (Line(
          points={{-148,-38},{108,-38},{108,-99.1}},
                                                  color={127,127,0}));
    connect(Cl_E.solution, blood_erythrocytes.solution) annotation (Line(
        points={{0,-38},{0,-99.1},{108,-99.1}},
        color={127,127,0}));
    connect(HCO3_E.solution, blood_erythrocytes.solution) annotation (Line(
          points={{-24,-38},{108,-38},{108,-99.1}},
                                                color={127,127,0}));
    connect(Aquapirin.port_b, H2O_E.port_a) annotation (Line(
        points={{-168,-10},{-168,-28},{-164,-28}},
        color={158,66,200},
        thickness=1));
    connect(Aquapirin.port_a, H2O.port_a) annotation (Line(
        points={{-168,10},{-168,30},{-166,30}},
        color={158,66,200},
        thickness=1));
    connect(Band3.port_a, HCO3.port_a) annotation (Line(
        points={{-6,10},{-6,30},{-8,30}},
        color={158,66,200},
        thickness=1));
    connect(Band3.port_b, HCO3_E.port_a) annotation (Line(
        points={{-6,-10},{-6,-28},{-8,-28}},
        color={158,66,200},
        thickness=1));
    connect(Band3_.port_b, Cl_E.port_a) annotation (Line(
        points={{18,-10},{18,-28},{16,-28}},
        color={158,66,200},
        thickness=1));
    connect(Band3_.port_a, Cl.port_a) annotation (Line(
        points={{18,10},{18,30},{16,30}},
        color={158,66,200},
        thickness=1));
  connect(HCO3.solution, blood_plasma.solution) annotation (Line(
      points={{-24,20},{108,20},{108,12.88}},
      color={127,127,0}));
    connect(blood_plasma.solution, permeableUncharged.solution) annotation (Line(
        points={{108,12.88},{108,20},{162,20}},
        color={127,127,0}));
    connect(blood_erythrocytes.solution, permeableUncharged_E.solution)
      annotation (Line(
        points={{108,-99.1},{108,-38},{160,-38}},
        color={127,127,0}));
    connect(blood_erythrocytes.solution,chargedImpermeable_E. solution)
      annotation (Line(
        points={{108,-99.1},{108,-38},{140,-38},{140,-62},{148,-62}},
        color={127,127,0}));
    connect(permeableUncharged.port_a, leak.port_a) annotation (Line(
        points={{146,30},{140,30},{140,10}},
        color={158,66,200},
        thickness=1));
    connect(permeableUncharged_E.port_a, leak.port_b) annotation (Line(
        points={{144,-28},{140,-28},{140,-10}},
        color={158,66,200},
        thickness=1));
    connect(MCT_.port_a, Lac.port_a) annotation (Line(
        points={{78,10},{78,30},{76,30}},
        color={158,66,200},
        thickness=1));
    connect(MCT_.port_b, Lac_E.port_a) annotation (Line(
        points={{78,-10},{78,-28},{76,-28}},
        color={158,66,200},
        thickness=1));
    connect(Lac.solution, blood_plasma.solution) annotation (Line(
        points={{60,20},{108,20},{108,12.88}},
        color={127,127,0}));
    connect(blood_erythrocytes.solution, Lac_E.solution) annotation (Line(
        points={{108,-99.1},{108,-38},{60,-38}},
        color={127,127,0}));
    connect(H_E.solution, blood_erythrocytes.solution) annotation (Line(
        points={{34,-38},{108,-38},{108,-99.1}},
        color={127,127,0}));
    connect(H_E.port_a, MCT.port_b) annotation (Line(
        points={{50,-28},{52,-28},{52,-10}},
        color={158,66,200},
        thickness=1));
    connect(MCT.port_a, H.port_a) annotation (Line(
        points={{52,10},{52,30},{50,30}},
        color={158,66,200},
        thickness=1));
    connect(blood_plasma.solution, H.solution) annotation (Line(
        points={{108,12.88},{108,20},{34,20}},
        color={127,127,0}));
    connect(CO2.port_a, freeCO2.port_a) annotation (Line(
        points={{-40,30},{-36,30},{-36,10}},
        color={158,66,200},
        thickness=1));
    connect(freeCO2.port_b, CO2_E.port_a) annotation (Line(
        points={{-36,-10},{-36,-28},{-38,-28}},
        color={158,66,200},
        thickness=1));
    connect(blood_plasma.solution, CO2.solution) annotation (Line(
        points={{108,12.88},{108,20},{-56,20}},
        color={127,127,0}));
    connect(CO2_E.solution, blood_erythrocytes.solution) annotation (Line(
        points={{-54,-38},{108,-38},{108,-99.1}},
        color={127,127,0}));
    connect(blood_plasma.solution, O2.solution) annotation (Line(
        points={{108,12.88},{108,20},{100,20}},
        color={127,127,0}));
    connect(O2_E.solution, blood_erythrocytes.solution) annotation (Line(
        points={{100,-38},{108,-38},{108,-99.1}},
        color={127,127,0}));
    connect(O2_E.port_a, freeO2.port_b) annotation (Line(
        points={{116,-28},{118,-28},{118,-10}},
        color={158,66,200},
        thickness=1));
    connect(freeO2.port_a, O2.port_a) annotation (Line(
        points={{118,10},{118,30},{116,30}},
        color={158,66,200},
        thickness=1));
    connect(H2O.solution, K.solution) annotation (Line(
        points={{-150,20},{-96,20}},
        color={158,66,200}));
    connect(H2O.solution, Na.solution) annotation (Line(
        points={{-150,20},{-122,20}},
        color={158,66,200}));
    connect(H2O_E.solution, Na_E.solution) annotation (Line(
        points={{-148,-38},{-122,-38}},
        color={158,66,200}));
    connect(H2O_E.solution, K_E.solution) annotation (Line(
        points={{-148,-38},{-108,-38}},
        color={158,66,200}));
    connect(H2O_E.solution, H2PO4_E.solution) annotation (Line(
        points={{-148,-38},{-80,-38}},
        color={158,66,200}));
    connect(ADP_E.solution, K_E.solution) annotation (Line(
        points={{-110,-62},{-110,-38},{-108,-38}},
        color={158,66,200}));
    connect(ATP_E.solution, Na_E.solution) annotation (Line(
        points={{-122,-62},{-122,-38}},
        color={127,127,0}));
    connect(H2O.solution, H2PO4.solution) annotation (Line(
        points={{-150,20},{-82,20}},
        color={127,127,0}));
    connect(HPO4_E.solution, H2PO4_E.solution) annotation (Line(
        points={{-80,-62},{-110,-62},{-110,-38},{-80,-38}},
        color={127,127,0}));
    connect(HPO4.solution, H2PO4.solution) annotation (Line(
        points={{-82,42},{-82,20}},
        color={158,66,200}));
    connect(albumin.solution, permeableUncharged.solution) annotation (Line(
        points={{112,76},{92,76},{92,20},{162,20}},
        color={127,127,0}));
    connect(globulins.solution, permeableUncharged.solution) annotation (Line(
        points={{146,76},{92,76},{92,20},{162,20}},
        color={127,127,0}));
    connect(Ca.solution, CO2.solution) annotation (Line(
        points={{-108,42},{-82,42},{-82,20},{-56,20}},
        color={127,127,0}));
    connect(Mg.solution, blood_erythrocytes.solution) annotation (Line(
        points={{-108,-84},{-108,-38},{108,-38},{108,-99.1}},
        color={127,127,0}));
    connect(DPG.solution, permeableUncharged_E.solution) annotation (Line(
        points={{124,-94},{140,-94},{140,-38},{160,-38}},
        color={158,66,200}));
    connect(hemoglobin.solution, permeableUncharged_E.solution) annotation (Line(
        points={{90,-94},{140,-94},{140,-38},{160,-38}},
        color={127,127,0}));
    connect(GSH.solution, permeableUncharged_E.solution) annotation (Line(
        points={{160,-94},{140,-94},{140,-38},{160,-38}},
        color={127,127,0}));
    annotation ( Documentation(info="<html>
<p>Blood eqiulibrium across erythrocyte membrane bewteen blood plasma and intracellular fluid of erythrocytes.</p>
<p>Data of blood status are from:</p>
<p>Raftos, J.E., Bulliman, B.T. and Kuchel, P.W. Evaluation of an electrochemical model of erythrocyte pH buffering using 31P nuclear magnetic resonance data. <i>The Journal of general physiology</i> 1990;95(6):1183-1204. </p>
</html>",  revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
      experiment(StopTime=1e-008),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},{180,100}}),
                         graphics));
  end RedCellMembrane;

  model LeadAcidBattery
    "The electrochemical cell: PbSO4(s) | Pb(s) | HSO4-(aq) , H+(aq) | PbO2(s) | PbSO4(s) + 2 H2O"
   extends Modelica.Icons.Example;

    Chemical.Components.Solution anode(ElectricGround=false)
      annotation (Placement(transformation(extent={{24,-76},{58,32}})));

    Chemical.Components.Solution cathode(ElectricGround=false)
      annotation (Placement(transformation(extent={{-80,-78},{-46,30}})));

    Chemical.Components.Solution solution1(ElectricGround=false)
      annotation (Placement(transformation(extent={{-26,-80},{2,20}})));

    Chemical.Components.Substance Pb(substanceData=Substances.Lead_solid,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{50,-66},{30,-46}})));
    Chemical.Components.Substance HSO4(substanceData=Substances.HydrogenSulfate_aqueous,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-2,-70},{-22,-50}})));
    Chemical.Components.Substance PbSO4_(substanceData=Substances.LeadSulfate_solid,
        amountOfSubstance_start=0.01)
      annotation (Placement(transformation(extent={{50,-32},{30,-12}})));
    Chemical.Components.Substance H(substanceData=Substances.Proton_aqueous,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-2,-42},{-22,-22}})));
    Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
      annotation (Placement(transformation(extent={{-32,72},{-12,92}})));
    Chemical.Components.Reaction electrodeReaction(
      nP=2,
      nS=4,
      s={1,1,3,2},
      p={1,2}) annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=90,
          origin={-36,-14})));
    Chemical.Components.Reaction electrodeReaction1(
      nS=2,
      nP=3,
      p={1,1,2}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={14,-16})));

  Chemical.Components.ElectronTransfer electrone
      annotation (Placement(transformation(extent={{50,2},{30,22}})));
  Chemical.Components.ElectronTransfer electrone1
      annotation (Placement(transformation(extent={{-72,-38},{-52,-18}})));
    Chemical.Components.Substance PbO2(substanceData=Substances.LeadDioxide_solid,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent=
              {{-10,-10},{10,10}}, origin={-60,-58})));
    Chemical.Components.Substance H2O(substanceData=Substances.Water_liquid,
        amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{-2,-8},{-22,12}})));
    Chemical.Components.Substance PbSO4(substanceData=Substances.LeadSulfate_solid,
        amountOfSubstance_start=0.01) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}}, origin={-60,6})));

  Modelica.Electrical.Analog.Basic.Ground ground
    annotation (Placement(transformation(extent={{16,30},{36,50}})));
  Modelica.Electrical.Analog.Basic.Resistor resistor(R=1)
    annotation (Placement(transformation(extent={{-14,40},{6,60}})));
  Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
    annotation (Placement(transformation(extent={{-56,40},{-36,60}})));
  equation
    connect(Pb.port_a, electrodeReaction1.substrates[1]) annotation (Line(
        points={{30,-56},{15.5,-56},{15.5,-26},{16,-26}},
        color={158,66,200},
        thickness=0.5));
    connect(HSO4.port_a, electrodeReaction1.substrates[2]) annotation (Line(
        points={{-22,-60},{12,-60},{12,-26}},
        color={158,66,200},
        thickness=0.5));
    connect(PbSO4_.port_a, electrodeReaction1.products[1]) annotation (Line(
        points={{30,-22},{26,-22},{26,-2},{16,-2},{16,-6},{16.6667,-6}},
        color={158,66,200},
        thickness=0.5));
    connect(HSO4.solution, solution1.solution) annotation (Line(
        points={{-6,-70},{-6,-70},{-6,-78},{-3.6,-78},{-3.6,-79}},
        color={127,127,0}));
    connect(H.solution, solution1.solution) annotation (Line(points={{-6,-42},
          {-6,-78},{-3.6,-78},{-3.6,-79}},
                                       color={127,127,0}));
    connect(H2O.solution, solution1.solution) annotation (Line(
        points={{-6,-8},{-6,-79},{-3.6,-79}},
        color={127,127,0}));
    connect(electrodeReaction.products[1], PbSO4.port_a) annotation (Line(
        points={{-38,-4},{-38,6},{-50,6}},
        color={158,66,200},
        thickness=0.5));
    connect(electrodeReaction.products[2], H2O.port_a) annotation (Line(
        points={{-34,-4},{-34,-4},{-34,2},{-22,2}},
        color={158,66,200},
        thickness=0.5));
    connect(PbO2.port_a, electrodeReaction.substrates[1]) annotation (Line(
        points={{-50,-58},{-36,-58},{-36,-24},{-39,-24}},
        color={158,66,200},
        thickness=0.5));
    connect(HSO4.port_a, electrodeReaction.substrates[2]) annotation (Line(
        points={{-22,-60},{-34,-60},{-34,-24},{-37,-24}},
        color={158,66,200},
        thickness=0.5));
    connect(H.port_a, electrodeReaction.substrates[3]) annotation (Line(
        points={{-22,-32},{-32,-32},{-32,-24},{-35,-24}},
        color={158,66,200},
        thickness=0.5));
    connect(electrone1.port_a, electrodeReaction.substrates[4]) annotation (Line(
        points={{-52,-28},{-38,-28},{-38,-24},{-33,-24}},
        color={158,66,200},
        thickness=0.5));
    connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(
        points={{-22,-32},{2,-32},{2,2},{12,2},{12,-6},{14,-6}},
        color={158,66,200},
        thickness=0.5));
    connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
        points={{30,12},{14,12},{14,-6},{11.3333,-6}},
        color={158,66,200},
        thickness=0.5));
  connect(Pb.solution, anode.solution) annotation (Line(
      points={{46,-66},{46,-74.92},{51.2,-74.92}},
      color={127,127,0}));
  connect(PbSO4_.solution, anode.solution) annotation (Line(
      points={{46,-32},{46,-74.92},{51.2,-74.92}},
      color={127,127,0}));
  connect(PbO2.solution, cathode.solution) annotation (Line(
      points={{-66,-68},{-66,-70},{-60,-70},{-60,-76.92},{-52.8,-76.92}},
      color={127,127,0}));
  connect(electrone1.pin, voltageSensor.p) annotation (Line(
      points={{-72,-28},{-82,-28},{-82,50},{-64,50},{-64,82},{-32,82}},
      color={0,0,255}));
  connect(electrone.pin, voltageSensor.n) annotation (Line(
      points={{50,12},{50,50},{26,50},{26,82},{-12,82},{-12,82}},
      color={0,0,255}));
  connect(electrone.solution, anode.solution) annotation (Line(
      points={{46,2},{46,-74.92},{51.2,-74.92}},
      color={127,127,0}));
  connect(electrone.pin, ground.p) annotation (Line(
      points={{50,12},{50,50},{26,50}},
      color={0,0,255}));
  connect(electrone1.pin, currentSensor.p) annotation (Line(
      points={{-72,-28},{-82,-28},{-82,50},{-56,50}},
      color={0,0,255}));
  connect(currentSensor.n, resistor.p) annotation (Line(
      points={{-36,50},{-14,50}},
      color={0,0,255}));
  connect(resistor.n, electrone.pin) annotation (Line(
      points={{6,50},{50,50},{50,12}},
      color={0,0,255}));
  connect(PbSO4.solution, cathode.solution) annotation (Line(
      points={{-66,-4},{-66,-70},{-60,-70},{-60,-76.92},{-52.8,-76.92}},
      color={127,127,0}));
  connect(electrone1.solution, cathode.solution) annotation (Line(
      points={{-68,-38},{-66,-38},{-66,-70},{-60,-70},{-60,-76.92},{-52.8,
            -76.92}},
      color={127,127,0}));

    annotation (
    experiment(StopTime=49500), Documentation(revisions=
                      "<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>The lead-acid electrochemical cells are characterized by two chemical reactions:</p>
<table width=100%>
<tr><th>PbO2 + HSO4- + 3 H+ +2 e- &harr; PbSO4 + 2 H2O</th><td>(1)</td></tr>
<tr><th>Pb + HSO4- &harr; PbSO4 + H+ + 2 e-</th><td>(2)</td></tr>
</table>
<p>The building of one cell of a lead-acid battery starts with the definition of three solutions: two for the lead elec-trodes and one for the liquid-acid solution (Figure 1A). This can be done by dragging and dropping the library class &lsquo;Components.Solution&rsquo; into the diagram. We called the first instance &ldquo;cathode&rdquo;, the second &ldquo;solution&rdquo; and the last &ldquo;anode&rdquo;. We set the parameter &lsquo;Electri-calGround&rsquo; as &ldquo;false&rdquo; for all of these solutions in order to attain the possibility of non-zero voltages. Now we can specify the chemical substances inside the chemical solutions. We drag and drop the library class &lsquo;Compo-nents.Substance&rsquo; into the &ldquo;solution&rdquo; as chemical sub-stances (Figure 1B). H2O(liquid), H+(aqueous) and HSO4-(aqueous) representing the liquid aqueous solu-tion of sulfuric acid. PbSO4(solid) and PbO2(solid) are placed in the &ldquo;cathode&rdquo;, representing the elements of the positive electrode. The substances Pb(solid) and aP-bSO4(solid) are placed into the &ldquo;anode&rdquo;, representing the elements of the negative electrode. All of these sub-stances must be given unique names (e.g., &ldquo;PbSO4&rdquo; for the cathode and &ldquo;aPbSO4&rdquo; for the anode), because the Modelica language does not support two instances with the same name in a single class.</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/LeadAcidBatterry1.png\"/></p>
<p>Figure 1) The building of one electro-chemical cell of a lead-acid battery in four steps: A) adding chemical solutions, B) adding chemical substances, C) adding electron transfers and D) adding chemical reactions.</p>
<p>As mentioned above, the appropriate substance data for all these substances must be selected as predefined parametric records, e.g., &lsquo;Exam-ples.Substances.Water_liquid&rsquo;, &lsquo;.Lead_solid&rsquo;, &lsquo;.Lead_dioxide_solid&rsquo;, &lsquo;.Lead_sulfate_solid&rsquo;, and so on. The last, very special substance to be included is an electron. This class is called &lsquo;Compo-nents.ElectronTransfer&rsquo; and it must be added in order for each electrode to transfer electron from the chemical reaction to the electric circuit (Figure 1C). Each of these substances must be connected to the appropriate solu-tion using a solution port situated in the bottom of the component&rsquo;s icons to indicate that they are all mixed in the solution. By having all these substances, it is possi-ble to implement the chemical reactions. Dragging and dropping the library class &lsquo;Components.Reaction&rsquo; for both chemical reactions, and setting their parameters as an appropriate number of reactants, products and stoi-chiometry, allows the connection of each substance with the reaction, as expressed in reaction (1) and reaction (2). This setting can be done using the parameter dialog of the cathode chemical reaction (1) as there are four types of substrates (nS=4) with stoichiometric coeffi-cients: one for the first and second reactant, three for the third reactant and two for the fourth reactant (s={1,1,3,2}). There are also two types of products (nP=2) with stoichiometry: one for PbSO4 and two for water (p={1,2}), following the chemical scheme of the first chemical reaction above. After setting the number of reactants and products, it is possible to connect the substances with reactions. Each instance of reaction has an array of connectors for substrates and an array of con-nectors for products; the user must be very careful to connect each element of these arrays in the same order as defined by stoichiometric coefficients. This means that, for example, the water must be connected in index 2 to products of the first chemical reaction, because we had already selected the order of products by setting the array of stoichiometric coefficients in reaction (1). The chemical reaction (2) must be set analogically as nS=2, nP=3, p={1,1,2} with connections of substance ports of Pb to substrate[1], HSO4- to substrate[2], PbSO4 to prod-uct[1], H+ to product[2] and e- to product[3], as repre-sented in Figure 1D.</p>
<p>The electrochemical cell has already been imple-mented at this stage. However, the simulation requires the initial state of substances, which for the fully charged battery means that almost all elements of the cathode are PbO2 and almost all elements of the anode are Pb. In this state, the sulfuric acid can be concen-trated, which increases the effectiveness of the electro-chemical cell. To set this state, it is possible to just dou-ble-click on PbO2 and Pb and set the amount, e.g., 1mol. To set the pure concentrated sulfuric acid we can also set the amount of SO4- and H+ as 1mol. This fully charged ideal state is ready to simulate when it is con-nected to the electric ground via one of the electric ports of the one electron transfer component.</p>
<p>These batteries can be connected to any electrical cir-cuit that is slowly discharging. For example, if we only connect the simple electric resistance of 1 Ohm as ex-pressed in Figure 1D, then the simulation of the dis-charging process over 13 hours and 45 minutes gives the results of electric current and electric potential, as can be seen in Figure 2. The exchange of the resistor with a voltage source can simulate the charging process for a discharged cell.</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/LeadAcidBatterry2.png\"/></p>
<p>Figure 2) Discharging simulation of the lead-acid battery cell from Figure 2D, with the initial amount of substances as described in the text.</p>
</html>"));
  end LeadAcidBattery;

  package AcidBase

    model WaterSelfIonization "H2O  <->  OH-   +   H+ "
      import Chemical;
        extends Modelica.Icons.Example;

      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-72,2},{76,96}})));
      Chemical.Components.Solution solution1
        annotation (Placement(transformation(extent={{-76,-98},{72,-4}})));
      Chemical.Components.Substance H3O(amountOfSubstance_start=1e-7,
          substanceData=Chemical.Examples.Substances.Hydronium_aqueous)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin=
               {30,70})));
      Chemical.Components.Substance OH(amountOfSubstance_start=1e-7,
          substanceData=Chemical.Examples.Substances.Hydroxide_aqueous)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin=
               {30,26})));
      Chemical.Components.Substance H2O(amountOfSubstance_start=1,
          substanceData=Chemical.Examples.Substances.Water_liquid) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={-30,46})));
      Chemical.Components.Reaction waterDissociation(nP=2, s={2})
        annotation (Placement(transformation(extent={{-12,36},{8,56}})));
            Real pH, pH_;
      Chemical.Components.Substance H_(amountOfSubstance_start=1e-7,
          substanceData=Chemical.Examples.Substances.Proton_aqueous)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin=
               {28,-30})));
      Chemical.Components.Substance OH_(amountOfSubstance_start=1e-7,
          substanceData=Chemical.Examples.Substances.Hydroxide_aqueous)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin=
               {28,-76})));
      Chemical.Components.Substance H2O_(amountOfSubstance_start=1,
          substanceData=Chemical.Examples.Substances.Water_liquid) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={-32,-56})));
      Chemical.Components.Reaction waterDissociation_(nP=2)
        annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));

    equation
      pH = -log10( H3O.a);

      pH_ = -log10( H_.a);

      connect(OH.port_a, waterDissociation.products[1]) annotation (Line(
          points={{20,26},{16,26},{16,44},{8,44}},
          color={158,66,200},
          thickness=1));
      connect(waterDissociation.products[2], H3O.port_a) annotation (Line(
          points={{8,48},{16,48},{16,70},{20,70}},
          color={158,66,200},
          thickness=1));
      connect(H2O.port_a, waterDissociation.substrates[1]) annotation (Line(
          points={{-20,46},{-12,46}},
          color={158,66,200},
          thickness=1));
      connect(OH_.port_a,waterDissociation_. products[1]) annotation (Line(
          points={{18,-76},{14,-76},{14,-58},{6,-58}},
          color={158,66,200},
          thickness=1));
      connect(waterDissociation_.products[2], H_.port_a) annotation (Line(
          points={{6,-54},{14,-54},{14,-30},{18,-30}},
          color={158,66,200},
          thickness=1));
      connect(H2O_.port_a,waterDissociation_. substrates[1]) annotation (Line(
          points={{-22,-56},{-14,-56}},
          color={158,66,200},
          thickness=1));
      connect(H2O.solution, solution.solution) annotation (Line(
          points={{-36,36},{46.4,36},{46.4,2.94}},
          color={127,127,0}));
      connect(OH.solution, solution.solution) annotation (Line(
          points={{36,16},{36,2.94},{46.4,2.94}},
          color={127,127,0}));
      connect(H3O.solution, solution.solution) annotation (Line(
          points={{36,60},{36,2.94},{46.4,2.94}},
          color={127,127,0}));
      connect(H2O_.solution, solution1.solution) annotation (Line(
          points={{-38,-66},{42.4,-66},{42.4,-97.06}},
          color={127,127,0}));
      connect(OH_.solution, solution1.solution) annotation (Line(
          points={{34,-86},{34,-97.06},{42.4,-97.06}},
          color={127,127,0}));
      connect(H_.solution, solution1.solution) annotation (Line(
          points={{34,-40},{34,-97.06},{42.4,-97.06}},
          color={127,127,0}));
      annotation ( Documentation(info="<html>
<p>Self-ionization of water.</p>
<p>Ions difference (SID) in water causes the acidity/basicity, where pH = -log10(aH+). An activity of hydrogen ions aH+ is approximated with concentration (mol/l) of the oxonium cations H3O+.</p>
<pre><b>plotExpression(apply(-log10(WaterSelfIonization.H3O.solute)),&nbsp;false,&nbsp;&quot;pH&quot;,&nbsp;1);</b></pre>
<p><br>The titration slope der(pH)/der(SID)=1.48e+6 1/(mol/L) at pH=7.4.</p>
</html>",    revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=1));
    end WaterSelfIonization;

    model CarbonDioxideInWater "CO2 as alone acid-base buffer"
      import Chemical;
        extends Modelica.Icons.Example;
      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,46}})));
      Chemical.Components.Substance HCO3(substanceData=Chemical.Examples.Substances.Bicarbonate_aqueous)
        annotation (Placement(transformation(extent={{-16,-4},{4,16}})));
      Chemical.Components.Reaction HendersonHasselbalch(nP=2, nS=2,
      useKineticsInput=false) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
        annotation (Placement(transformation(extent={{-48,-6},{-28,14}})));
    Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
      substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
      PartialPressure=5332.8954966,
      TotalPressure=101325.0144354) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-60,86})));
      Chemical.Components.Substance H(substanceData=Chemical.Examples.Substances.Proton_aqueous,
          amountOfSubstance_start=3e-8) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}}, origin={-6,-38})));
      Chemical.Components.GasSolubility gasSolubility
        annotation (Placement(transformation(extent={{-70,36},{-50,56}})));
                                            /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
      Chemical.Components.Substance CO2_liquid(substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous)
        annotation (Placement(transformation(extent={{-82,-6},{-62,14}})));
      Chemical.Components.Substance CO3(substanceData=Chemical.Examples.Substances.Carbonate_aqueous)
        annotation (Placement(transformation(extent={{70,-2},{50,18}})));
      Chemical.Components.Reaction c2(nP=2, nS=1)
        "K=10^(-10.33 + 3), dH=14.9kJ/mol"
        annotation (Placement(transformation(extent={{16,-4},{36,16}})));
      Chemical.Components.Substance H2O(  substanceData=
            Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=55.507)
        annotation (Placement(transformation(extent={{-82,-50},{-62,-30}})));
      Real pH;

    equation
      pH = -log10( H.a);

      connect(CO2_gas.port_a, gasSolubility.gas_port) annotation (Line(
          points={{-60,76},{-60,56}},
          color={158,66,200},
          thickness=1));
      connect(gasSolubility.liquid_port, CO2_liquid.port_a) annotation (Line(
          points={{-60,36},{-60,4},{-62,4}},
          color={158,66,200},
          thickness=1));
      connect(HendersonHasselbalch.products[1], H.port_a) annotation (Line(
          points={{-28,2},{-22,2},{-22,-38},{4,-38}},
          color={158,66,200},
          thickness=1));
      connect(HendersonHasselbalch.products[2], HCO3.port_a) annotation (Line(
          points={{-28,6},{4,6}},
          color={158,66,200},
          thickness=1));
      connect(HCO3.port_a, c2.substrates[1]) annotation (Line(
          points={{4,6},{16,6}},
          color={158,66,200},
          thickness=1));
      connect(c2.products[1], H.port_a) annotation (Line(
          points={{36,4},{44,4},{44,-38},{4,-38}},
          color={158,66,200},
          thickness=1));
      connect(c2.products[2], CO3.port_a) annotation (Line(
          points={{36,8},{48,8},{48,8},{50,8}},
          color={158,66,200},
          thickness=1));
      connect(CO2_liquid.port_a, HendersonHasselbalch.substrates[2]) annotation (
          Line(
          points={{-62,4},{-62,6},{-48,6}},
          color={158,66,200},
          thickness=1));
      connect(H2O.port_a, HendersonHasselbalch.substrates[1]) annotation (Line(
          points={{-62,-40},{-56,-40},{-56,2},{-48,2}},
          color={158,66,200},
          thickness=1));
      connect(CO2_liquid.solution, solution.solution) annotation (Line(
          points={{-78,-6},{-78,-98.54},{60,-98.54}},
          color={127,127,0}));
      connect(H2O.solution, solution.solution) annotation (Line(points={{-78,-50},
            {-78,-98.54},{60,-98.54}},color={127,127,0}));
      connect(HCO3.solution, solution.solution) annotation (Line(points={{-12,-4},
            {-12,-98.54},{60,-98.54}},color={127,127,0}));
      connect(H.solution, solution.solution) annotation (Line(points={{-12,-48},
            {-12,-98.54},{60,-98.54}},color={127,127,0}));
      connect(CO3.solution, solution.solution) annotation (Line(points={{66,-2},
            {66,-98.54},{60,-98.54}}, color={127,127,0}));
      annotation ( Documentation(info="<html>
<p>CO2 solution in water without any other acid-base buffers.</p>
<pre><b>plotExpression(apply(-log10(CarbonDioxideInWater.H3O.solute)),&nbsp;false,&nbsp;&quot;pH&quot;,&nbsp;1);</b></pre>
<p><br>Please note, that OH- (and CO3^-2) can be neglected from electroneutrality calculation, because of very small concentrations (in physiological pH) anyway. </p>
<p>And if SID&gt;0 then also H3O+ can be also neglected from electroneutrality, because only bicarbonate anions HCO3- (or CO3^-2) are needed there to balance the electroneutrality.</p>
<p><br>The partial pressure of CO2 in gas are input parameter. Outputs are an amount of free dissolved CO2 in liquid and an amount of HCO3-.</p>
<p><br>The titration slope der(pH)/der(SID)=17.5 1/(mol/L) at pH=7.4 and pCO2=40 mmHg.</p>
<p><br>Molar heat of formation (aqueous):</p>
<p>CO2:        -413.5 kJ/mol  (gas: -393.5 kJ/mol )</p>
<p>H2O:        -285.8 kJ/mol</p>
<p>HCO3-:        -692.0 kJ/mol</p>
<p>CO3^-2:        -677.1 kJ/mol</p>
<p><br>Enthalphy of reaction H2O + CO2 &lt;-&gt; HCO3- + H+  :         7.3 kJ/mol</p>
<p>Enthalphy of reaction HCO3- &lt;-&gt; CO3^-2 + H+  :        14.9 kJ/mol</p>
</html>",    revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.02));
    end CarbonDioxideInWater;

    model Phosphate
      import Chemical;
        extends Modelica.Icons.Example;

      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-98,-100},{100,100}})));

       Chemical.Components.Substance H(amountOfSubstance_start=55.6*10^(-7.4),
          substanceData=Chemical.Examples.Substances.Proton_aqueous)
        "hydrogen ions activity" annotation (Placement(transformation(extent={{
                -10,-10},{10,10}}, origin={28,-14})));

      Chemical.Components.Substance H3PO4(amountOfSubstance_start=1e-08,
          substanceData=Chemical.Examples.Substances.PhosphoricAcid_aqueous)
        annotation (Placement(transformation(extent={{-90,-58},{-70,-38}})));
      Chemical.Components.Substance H2PO4(amountOfSubstance_start=0.0005,
          substanceData=Chemical.Examples.Substances.DihydrogenPhosphate_aqueous)
        annotation (Placement(transformation(extent={{-40,-58},{-20,-38}})));
      Chemical.Components.Substance HPO4(substanceData=Chemical.Examples.Substances.HydrogenPhosphate_aqueous,
          amountOfSubstance_start=0.0006)
        annotation (Placement(transformation(extent={{16,-58},{36,-38}})));
      Chemical.Components.Substance PO4(substanceData=Chemical.Examples.Substances.Phosphate_aqueous,
          amountOfSubstance_start=1e-08)
        annotation (Placement(transformation(extent={{92,-58},{72,-38}})));

      Chemical.Components.Reaction chemicalReaction(nP=2) "10^(-1.915 + 3)"
        annotation (Placement(transformation(extent={{-66,-58},{-46,-38}})));
      Chemical.Components.Reaction chemicalReaction1(nP=2) "10^(-6.66 + 3)"
        annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
      Chemical.Components.Reaction chemicalReaction2(nP=2) "10^(-11.78 + 3)"
        annotation (Placement(transformation(extent={{44,-58},{64,-38}})));

      Chemical.Components.Substance
                           H2O(      substanceData=
            Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=55.508)
                                   annotation (Placement(transformation(extent={{-10,
                -10},{10,10}}, origin={58,-76})));
    equation
      connect(H3PO4.port_a, chemicalReaction.substrates[1]) annotation (Line(
          points={{-70,-48},{-66,-48}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction.products[1], H2PO4.port_a) annotation (Line(
          points={{-46,-50},{-42,-50},{-42,-48},{-20,-48}},
          color={107,45,134},
          thickness=1));
      connect(H2PO4.port_a, chemicalReaction1.substrates[1]) annotation (Line(
          points={{-20,-48},{-14,-48}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction1.products[1], HPO4.port_a) annotation (Line(
          points={{6,-50},{16,-50},{16,-48},{36,-48}},
          color={107,45,134},
          thickness=1));
      connect(HPO4.port_a, chemicalReaction2.substrates[1]) annotation (Line(
          points={{36,-48},{44,-48}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction2.products[1], PO4.port_a) annotation (Line(
          points={{64,-50},{74,-50},{74,-48},{72,-48}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction.products[2], H.port_a) annotation (Line(
          points={{-46,-46},{-44,-46},{-44,-32},{38,-32},{38,-14}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction1.products[2], H.port_a) annotation (Line(
          points={{6,-46},{14,-46},{14,-32},{38,-32},{38,-14}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction2.products[2], H.port_a) annotation (Line(
          points={{64,-46},{66,-46},{66,-32},{38,-32},{38,-14}},
          color={107,45,134},
          thickness=1));
      connect(H3PO4.solution, solution.solution) annotation (Line(
          points={{-86,-58},{-46,-58},{-46,-98},{60.4,-98}}));
      connect(H2PO4.solution, solution.solution) annotation (Line(points={{-36,-58},
            {-36,-88},{60.4,-88},{60.4,-98}}));
      connect(HPO4.solution, solution.solution) annotation (Line(points={{20,-58},
            {22,-58},{22,-88},{60.4,-88},{60.4,-98}}));
      connect(PO4.solution, solution.solution) annotation (Line(points={{88,-58},
            {88,-88},{60.4,-88},{60.4,-98}}));
      connect(H.solution, solution.solution) annotation (Line(points={{22,-24},
            {22,-88},{60.4,-88},{60.4,-98}}));
    connect(chemicalReaction.substrates[1], H3PO4.port_a) annotation (Line(
        points={{-66,-48},{-70,-48}},
        color={158,66,200},
        thickness=1));
    connect(H2O.solution, solution.solution) annotation (Line(
        points={{52,-86},{52,-98},{60.4,-98}},
        color={158,66,200}));
      annotation ( Documentation(info="<html>
</html>",    revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.05));
    end Phosphate;

    model AlbuminTitration "Figge-Fencl model (22. Dec. 2007)"
      extends Modelica.Icons.Example;

      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-104,-100},{96,100}})));

    Chemical.Sources.Buffer H(substanceData=Substances.Proton_aqueous, a_start=
            10^(-7.4)) "hydrogen ions activity" annotation (Placement(
            transformation(extent={{10,-10},{-10,10}}, origin={14,36})));

      constant Integer n=218 "Number of weak acid group in albumin molecule";
      constant Real pKAs[n]=cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
        "acid dissociation constants";
      constant Real K[n]=fill(10.0, n) .^ (-pKAs);
      constant Real DfG[n]= Modelica.Constants.R*(298.15)*log(K);

      Chemical.Components.Substance A[n](
        each amountOfSubstance_start=0.00033, substanceData(each z=-1))
        "deprotonated acid groups"
        annotation (Placement(transformation(extent={{26,-16},{6,4}})));
      Chemical.Components.Reaction react[n](each nP=2)
        annotation (Placement(transformation(extent={{-44,-2},{-24,18}})));

      Chemical.Components.Substance HA[n](substanceData(DfG_25degC_1bar=DfG), each amountOfSubstance_start=
           0.00033) "protonated acid groups"
        annotation (Placement(transformation(extent={{-78,-2},{-58,18}})));

      Chemical.Components.Substance H2O(substanceData=Substances.Water_liquid,
          amountOfSubstance_start=55.508) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}}, origin={62,-68})));
    equation
      connect(react.products[1], A.port_a) annotation (Line(
          points={{-24,6},{-12,6},{-12,-6},{6,-6}},
          color={107,45,134},
          thickness=1));
      for i in 1:n loop
        connect(react[i].products[2], H.port_a) annotation (Line(
            points={{-24,10},{-14,10},{-14,36},{4,36}},
            color={107,45,134},
            thickness=1));
        connect(HA[i].solution, solution.solution) annotation (Line(
          points={{-74,-2},{-74,-86},{56,-86},{56,-98}},
          color={127,127,0}));
        connect(A[i].solution, solution.solution) annotation (Line(
          points={{22,-16},{22,-86},{56,-86},{56,-98}},
          color={127,127,0}));
      end for;
      connect(HA.port_a, react.substrates[1]) annotation (Line(
          points={{-58,8},{-44,8}},
          color={107,45,134},
          thickness=1));

    connect(solution.solution, H2O.solution) annotation (Line(
        points={{56,-98},{56,-78}},
        color={127,127,0}));
    connect(H.solution, solution.solution) annotation (Line(
        points={{20,26},{20,14},{36,14},{36,-98},{56,-98}},
        color={127,127,0}));
      annotation ( Documentation(revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>The titration slope der(pH)/der(SID)=185 1/(mol/L) at pH=7.4 and tAlb=0.66 mmol/l.</p>
<p>Data and model is described in</p>
<p><font style=\"color: #222222; \">Jame Figge: Role of non-volatile weak acids (albumin, phosphate and citrate). In: Stewart&apos;s Textbook of Acid-Base, 2nd Edition, John A. Kellum, Paul WG Elbers editors, &nbsp;AcidBase org, 2009, pp. 216-232.</font></p>
</html>"),
        experiment(StopTime=1.6));
    end AlbuminTitration;

    model CarbonDioxideInBlood
      import Chemical;
        extends Modelica.Icons.Example;

      parameter Real KC=10;//e-6 "Slow down factor";

      Chemical.Components.Solution blood_erythrocytes(ElectricGround=false,
          temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-98},{100,-38}})));
      Chemical.Components.Solution blood_plasma(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,4},{100,56}})));

      Chemical.Components.Substance HCO3(
          substanceData=Chemical.Examples.Substances.Bicarbonate_blood,
          amountOfSubstance_start=0.024)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin=
               {18,24})));
    Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
      substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
      TotalPressure(displayUnit="mmHg") = 101325.0144354,
      PartialPressure(displayUnit="mmHg") = 5332.8954966,
      usePartialPressureInput=true,
        Temperature=310.15)
                          annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-84,84})));
      Chemical.Components.GasSolubility gasSolubility(KC=KC)
        annotation (Placement(transformation(extent={{-94,48},{-74,68}})));

      Chemical.Components.Substance CO2(substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous,
          amountOfSubstance_start=0.00148) "Free dissolved CO2 in plasma"
        annotation (Placement(transformation(extent={{-88,28},{-68,48}})));
      Chemical.Components.Substance H2O(substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=51.6159)
        annotation (Placement(transformation(extent={{-60,12},{-40,32}})));
      Chemical.Components.Substance HCO3_E(
          substanceData=Chemical.Examples.Substances.Bicarbonate_blood,
          amountOfSubstance_start=0.0116)
        annotation (Placement(transformation(extent={{28,-60},{8,-40}})));
      Chemical.Components.Reaction HendersonHasselbalch1(nP=2, nS=2,
      KC=KC) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
        annotation (Placement(transformation(extent={{-26,-68},{-6,-48}})));
      Chemical.Components.Substance CO2_E(substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous,
          amountOfSubstance_start=0.0011) "Free dissolved CO2 in erythrocyte"
        annotation (Placement(transformation(extent={{-90,-82},{-70,-62}})));
      Chemical.Components.Substance H2O_E(substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=38.4008)
        annotation (Placement(transformation(extent={{-60,-62},{-40,-42}})));
      Chemical.Components.Substance Cl_E(
          substanceData=Chemical.Examples.Substances.Chloride_aqueous,
          amountOfSubstance_start=0.0499)
        annotation (Placement(transformation(extent={{68,-60},{48,-40}})));
      Chemical.Components.Substance Cl(
          substanceData=Chemical.Examples.Substances.Chloride_aqueous,
          amountOfSubstance_start=0.103)
        annotation (Placement(transformation(extent={{68,20},{48,40}})));

      Real pH_e, pH_p;

      Chemical.Components.Membrane aquaporin(KC=KC) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-34,-16})));
      Chemical.Components.Membrane Band3_HCO3(KC=KC) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={4,-16})));
      Chemical.Components.Membrane Band3_Cl(useKineticsInput=false, KC=KC)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={46,-16})));
    Chemical.Sources.Buffer H_E(
      substanceData=Chemical.Examples.Substances.Proton_aqueous,
      BufferValue=0.063,
        a_start=10^(-7.2))
      annotation (Placement(transformation(extent={{48,-84},{30,-66}})));
      Modelica.Blocks.Sources.Clock clock(offset=5000)
        annotation (Placement(transformation(extent={{-54,62},{-34,82}})));
      Chemical.Components.Substance others_E(      substanceData(
          density=(1.045 - 0.695523)*1000/(1 - 0.697583),
          References={"erythrocyte intracellular fluid density 1045kg/m3"},
          MolarWeight=(1.045 - 0.695523)/(38.7*(1 - 0.994648) - 0.0499 - 0.0116
               - 0.00123)), amountOfSubstance_start=0.1444)
        annotation (Placement(transformation(extent={{68,-88},{88,-68}})));
      Chemical.Components.Substance others_P(   substanceData(
          References={
              "to reach plasma density 1024 kg/m3 and plasma volume 1 liter"},
          density=(1.024 - 0.933373)*1000/(1 - 0.936137),
          MolarWeight=(1.024 - 0.933373)/(51.8*(1 - 0.994648) - 0.103 - 0.024
               - 0.0017)), amountOfSubstance_start=0.1487)
        annotation (Placement(transformation(extent={{70,14},{90,34}})));
      Chemical.Components.Diffusion diffusion annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-66,-16})));
    Chemical.Sources.Buffer H(
      substanceData=Chemical.Examples.Substances.Proton_aqueous,
      BufferValue=0.0077,
        a_start=10^(-7.4))
        "buffer value 7.7 mmol/L for plasma is from (O. Siggaard-Andersen 1995)"
      annotation (Placement(transformation(extent={{38,38},{20,56}})));
      Chemical.Components.Reaction HendersonHasselbalch2(nP=2, nS=2,
        KC=(1e-10)*KC) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
        annotation (Placement(transformation(extent={{-26,26},{-6,46}})));
    equation
      pH_p = -log10(H.a);
      pH_e = -log10(H_E.a);
      connect(HendersonHasselbalch1.products[1], HCO3_E.port_a) annotation (Line(
          points={{-6,-60},{2,-60},{2,-50},{8,-50}},
          color={107,45,134},
          thickness=0.5));
    connect(CO2_E.port_a, HendersonHasselbalch1.substrates[1]) annotation (
        Line(
        points={{-70,-72},{-36,-72},{-36,-60},{-26,-60}},
        color={107,45,134},
        thickness=0.5));
      connect(H2O_E.port_a, HendersonHasselbalch1.substrates[2]) annotation (Line(
          points={{-40,-52},{-34,-52},{-34,-56},{-26,-56}},
          color={158,66,200},
          thickness=0.5));
    connect(CO2.solution, blood_plasma.solution) annotation (Line(
        points={{-84,28},{-84,12},{60,12},{60,4.52}},
        color={127,127,0}));
    connect(H2O.solution, blood_plasma.solution)
      annotation (Line(points={{-56,12},{-56,12},{60,12},{60,10},{60,4},{60,
              4.52}},                                   color={127,127,0}));
    connect(Cl.solution, blood_plasma.solution) annotation (Line(
        points={{64,20},{64,12},{60,12},{60,4.52}},
        color={127,127,0}));
    connect(CO2_E.solution, blood_erythrocytes.solution) annotation (Line(
        points={{-86,-82},{-86,-88},{60,-88},{60,-97.4}},
        color={127,127,0}));
      connect(H2O_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{-56,-62},{-56,-88},{60,-88},{60,-97.4}},
                                                    color={127,127,0}));
      connect(Cl_E.solution, blood_erythrocytes.solution) annotation (Line(
          points={{64,-60},{64,-78},{60,-78},{60,-97.4}},
          color={127,127,0}));
      connect(HCO3_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{24,-60},{24,-88},{60,-88},{60,-97.4}},
                                                  color={127,127,0}));
    connect(gasSolubility.liquid_port, CO2.port_a) annotation (Line(
        points={{-84,48},{-84,38},{-68,38}},
        color={158,66,200},
        thickness=0.5));
      connect(aquaporin.port_b, H2O_E.port_a) annotation (Line(
          points={{-34,-26},{-34,-52},{-40,-52}},
          color={158,66,200},
          thickness=0.5));
    connect(aquaporin.port_a, H2O.port_a) annotation (Line(
        points={{-34,-6},{-34,22},{-40,22}},
        color={158,66,200},
        thickness=0.5));
      connect(Band3_HCO3.port_a, HCO3.port_a) annotation (Line(
          points={{4,-6},{4,24},{8,24}},
          color={158,66,200},
          thickness=0.5));
      connect(Band3_HCO3.port_b, HCO3_E.port_a) annotation (Line(
          points={{4,-26},{4,-50},{8,-50}},
          color={158,66,200},
          thickness=0.5));
      connect(Band3_Cl.port_b, Cl_E.port_a) annotation (Line(
          points={{46,-26},{46,-38},{46,-50},{48,-50}},
          color={158,66,200},
          thickness=0.5));
      connect(Band3_Cl.port_a, Cl.port_a) annotation (Line(
          points={{46,-6},{46,12},{46,30},{48,30}},
          color={158,66,200},
          thickness=0.5));
      connect(gasSolubility.gas_port, CO2_gas.port_a) annotation (Line(
          points={{-84,68},{-84,74}},
          color={158,66,200},
          thickness=0.5));
    connect(HCO3.solution, blood_plasma.solution) annotation (Line(
        points={{24,14},{24,12},{60,12},{60,8},{60,4},{60,4.52}},
        color={127,127,0}));
    connect(H_E.port_a, HendersonHasselbalch1.products[2]) annotation (Line(
        points={{30,-75},{4,-75},{4,-56},{-6,-56}},
        color={158,66,200},
        thickness=0.5));
    connect(blood_erythrocytes.solution, others_E.solution) annotation (Line(
        points={{60,-97.4},{60,-88},{72,-88}},
        color={127,127,0}));
    connect(blood_plasma.solution, others_P.solution) annotation (Line(
        points={{60,4.52},{60,4},{60,8},{60,12},{74,12},{74,14}},
        color={127,127,0}));
    connect(clock.y, CO2_gas.partialPressure) annotation (Line(
        points={{-33,72},{-24,72},{-24,98},{-84,98},{-84,94}},
        color={0,0,127}));
    connect(H_E.solution, blood_erythrocytes.solution) annotation (Line(
        points={{44.4,-84},{44,-84},{44,-88},{60,-88},{60,-97.4}},
        color={127,127,0}));
      connect(CO2_E.port_a, diffusion.port_b) annotation (Line(
          points={{-70,-72},{-66,-72},{-66,-26}},
          color={158,66,200},
          thickness=0.5));
      connect(CO2.port_a, diffusion.port_a) annotation (Line(
          points={{-68,38},{-66,38},{-66,-6}},
          color={158,66,200},
          thickness=0.5));
      connect(blood_plasma.solution, H.solution) annotation (Line(
          points={{60,4.52},{60,4.52},{60,12},{34,12},{34,38},{34.4,38}},
          color={127,127,0}));
      connect(CO2.port_a, HendersonHasselbalch2.substrates[2]) annotation (Line(
          points={{-68,38},{-26,38}},
          color={158,66,200},
          thickness=0.5));
      connect(H2O.port_a, HendersonHasselbalch2.substrates[1]) annotation (Line(
          points={{-40,22},{-34,22},{-34,34},{-26,34}},
          color={158,66,200},
          thickness=0.5));
      connect(HendersonHasselbalch2.products[1], HCO3.port_a) annotation (Line(
          points={{-6,34},{2,34},{2,24},{8,24}},
          color={158,66,200},
          thickness=0.5));
      connect(HendersonHasselbalch2.products[2], H.port_a) annotation (Line(
          points={{-6,38},{2,38},{2,48},{20,48},{20,47}},
          color={158,66,200},
          thickness=0.5));
      annotation ( Documentation(info="<html>
<p>The mature red blood cell (erythrocyte) is the simplest cell in the human body. Its primary function is the transportation of blood gases, such as oxygen O<sub>2</sub> (from the lungs to tissues) and carbon dioxide CO<sub>2</sub> (from tissues to the lungs). The chemical processes behind the gases&rsquo; transportation are complex because the capacity of water to transport their freely dissolved forms is very low. To transport sufficient amounts of O<sub>2</sub> and CO<sub>2</sub>, the gases must be chemically bound to hemoglobin such as described in (Matej&aacute;k, et al., 2015) and/or transported as different substances, which can be present in water in much higher concentrations than their freely dissolved forms allow. Therefore, to transport a sufficient amount of CO<sub>2</sub>, it must be changed to HCO<sub>3</sub><sup>-</sup> using the chemical reaction: </p>
<table width=100%><tr>
<th><p align=\"center\"><b>CO<sub>2</sub> + H<sub>2</sub>O &lt;-&gt; HCO<sub>3</sub><sup>-</sup> + H<sup>+</sup></b></p></th>
<td><p>(1)</p></td>
</tr>
</table>
<p><br>This reaction takes place mainly inside the red blood cell, because only here it is presented with the enzyme carbonic anhydrase. Therefore, the increase of total carbon dioxide content of blood in tissues and its decrease in lungs are always connected with the chloride shift between blood plasma and the intracellular fluid of erythrocytes, as represented in followin Figure: </p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/CO2inBlood.png\"/></p>
<p>Figure: Chloride shift with carbon dioxide hydration with assumption of non-bicarbonate linear acid-base buffering properties of plasma and erythrocytes. </p>
<p><br>The blood plasma and intracellular fluid are divided by the cellular membrane composed of a special, very compact lipid double-layer. A lipophobic compound (not soluble in lipids) cannot cross the membrane without special proteins called membrane channels. Even water molecules must have membrane channels (called aquaporins) in order to cross the cellular membrane. In addition, the chloride shift (also known as the Hamburger shift) is exchanging an aqueous chloride Cl<sup>-</sup> for an aqueous bicarbonate HCO<sub>3</sub><sup>-</sup> in both directions across the cellular membranes of red blood cells using the membrane channel &ldquo;Band 3&rdquo;. Each passive membrane channel only allows the equilibration of the electrochemical potentials of the specific permeable ions on both sides of membrane. The different electric potentials on each side of membrane allow their different concentrations to achieve equilibrium. </p>
<p>Conversely, the solution&rsquo;s equilibrium of different ions&rsquo; compositions on both sides of the membrane creates the measurable electric membrane potential. This process is not so intuitive, because even though neither solution needs to have an electric charge, there can be a non-zero electric potential for permeable ions. This potential for permeable ions at equilibrium is called the Nernst membrane potential and, in the Chemical library, it is a direct mathematical result of the equality of the electrochemical potential of the ion in both solutions. </p>
<p>The intracellular solution must be set at the possible nonzero electric potential (ElectricalGround=false) because, as a result, the membrane potential of the erythrocytes is calculated as -12mV, which agrees with experimental data by Gedde and Huestis (Gedde and Huestis, 1997) in the electrolytes&rsquo; setting by Raftos et al. (Raftos, et al., 1990). </p>
<p>In this way, it is possible to model more complex processes of a membrane where chemical reactions of active membrane channels or membrane receptors can both be used.&nbsp; </p>
<p><br>CO2 in blood with linear H+ non-bicarbonates buffering without binding to hemoglobin.</p>
<p>The buffer values 0.063 mmol/L commes from Siggaard-Andersen.</p>
</html>",    revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(
        StopTime=1000),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end CarbonDioxideInBlood;

    model AcidBaseBufferTest
        extends Modelica.Icons.Example;

      Chemical.Sources.Buffer buffer(
        substanceData(z=1.045),
        a_start=10^(-7.2),
        BufferValue=3)
        annotation (Placement(transformation(extent={{-50,4},{-30,24}})));
      Chemical.Components.Solution simpleSolution
        annotation (Placement(transformation(extent={{-104,-100},{96,100}})));
      Chemical.Sources.ExternalMoleFraction externalMoleFraction(substanceData=
            Substances.Proton_aqueous,                   MoleFraction=10^(-7.1))
        annotation (Placement(transformation(extent={{0,-46},{20,-26}})));
      Chemical.Components.Substance substance(substanceData=Substances.Water_liquid,
          amountOfSubstance_start=1)
        annotation (Placement(transformation(extent={{52,-82},{72,-62}})));
    equation
      connect(buffer.solution, simpleSolution.solution) annotation (Line(
          points={{-46,4},{-26,4},{-26,-98},{56,-98}},
          color={127,127,0}));
      connect(externalMoleFraction.port_a, buffer.port_a) annotation (Line(
          points={{20,-36},{40,-36},{40,10},{-30,10},{-30,14}},
          color={158,66,200}));
      connect(simpleSolution.solution, substance.solution) annotation (Line(
          points={{56,-98},{26,-98},{26,-82},{56,-82}},
          color={127,127,0}));
      annotation (                experiment(StopTime=0.05));
    end AcidBaseBufferTest;

    package Dev
      model RedCellMembrane
       // import Chemical;
        extends Modelica.Icons.Example;

        parameter Real KC=1;//e-6 "Slow down factor";

        Chemical.Components.Solution blood_erythrocytes(ElectricGround=false)
          annotation (Placement(transformation(extent={{-180,-100},{180,-10}})));
        Chemical.Components.Solution blood_plasma
          annotation (Placement(transformation(extent={{-180,12},{180,100}})));

        Chemical.Components.Substance HCO3(substanceData=Substances.Bicarbonate_blood,
            amountOfSubstance_start(displayUnit="mmol") = 0.024) annotation (
            Placement(transformation(extent={{-10,-10},{10,10}}, origin={-18,30})));

        Chemical.Components.Substance H2O(substanceData=Substances.Water_liquid,
            amountOfSubstance_start=51.8*0.994648)
          annotation (Placement(transformation(extent={{-146,44},{-166,64}})));
        Chemical.Components.Substance HCO3_E(substanceData=Substances.Bicarbonate_blood,
            amountOfSubstance_start(displayUnit="mmol") = 0.0116)
          annotation (Placement(transformation(extent={{-28,-38},{-8,-18}})));
        Chemical.Components.Substance H2O_E(substanceData=Substances.Water_liquid,
            amountOfSubstance_start=38.7*0.994648)
          annotation (Placement(transformation(extent={{-144,-38},{-164,-18}})));
        Chemical.Components.Substance Cl_E(substanceData=Substances.Chloride_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.0499)
          annotation (Placement(transformation(extent={{-4,-38},{16,-18}})));
        Chemical.Components.Substance Cl(substanceData=Substances.Chloride_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.103)
          annotation (Placement(transformation(extent={{-4,20},{16,40}})));
        Chemical.Components.Substance albumin(substanceData(
                MolarWeight=66.463,
                z=-17,
                density=1080), amountOfSubstance_start(displayUnit="mmol") = 0.0007)
              annotation (Placement(transformation(extent={{112,76},{92,96}})));
      //  Real pH_e; //,pH_p;

        Chemical.Components.Membrane Aquapirin(KC=KC) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-168,0})));
        Chemical.Components.Membrane Band3(KC=KC) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-6,0})));
        Chemical.Components.Membrane Band3_(useKineticsInput=false, KC=KC)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={18,0})));
        Chemical.Components.Substance permeableUncharged(amountOfSubstance_start(
              displayUnit="mmol") = 0.0118)
          annotation (Placement(transformation(extent={{166,20},{146,40}})));
        Chemical.Components.Substance permeableUncharged_E(amountOfSubstance_start(
              displayUnit="mmol") = 0.00903, substanceData(MolarWeight=0.1))
          annotation (Placement(transformation(extent={{164,-38},{144,-18}})));
        Chemical.Components.Substance chargedImpermeable_E(amountOfSubstance_start(
              displayUnit="mmol") = 0.0165, substanceData(MolarWeight=1))
          annotation (Placement(transformation(extent={{144,-62},{164,-42}})));
        Chemical.Components.Membrane leak(useKineticsInput=false, KC=KC)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={140,0})));
        Chemical.Components.Substance Lac_E(substanceData=Substances.Chloride_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.00062)
          annotation (Placement(transformation(extent={{56,-38},{76,-18}})));
        Chemical.Components.Substance Lac(substanceData=Substances.Chloride_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.00131)
          annotation (Placement(transformation(extent={{56,20},{76,40}})));
        Chemical.Components.Membrane MCT_(useKineticsInput=false, KC=KC)
          "Monocarboxylate transporters" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={78,0})));
        Chemical.Components.Substance H_E(substanceData=Substances.Proton_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 38.7*10^(-7.2)) "H+"
          annotation (Placement(transformation(extent={{30,-38},{50,-18}})));
        Chemical.Components.Substance H(substanceData=Substances.Proton_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 51.8*10^(-7.4))
          "H+ in plasma"
          annotation (Placement(transformation(extent={{30,20},{50,40}})));
        Chemical.Components.Membrane MCT(useKineticsInput=false, KC=KC)
          "Monocarboxylate transporters" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={52,0})));
        Chemical.Components.Substance CO2(substanceData=Substances.CarbonDioxide_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.00167)
          "free dissolved unbound CO2"
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        Chemical.Components.Substance CO2_E(substanceData=Substances.CarbonDioxide_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.00125)
          "free dissolved unbound CO2"
          annotation (Placement(transformation(extent={{-58,-38},{-38,-18}})));
        Chemical.Components.Membrane freeCO2(KC=KC) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-38,2})));
        Chemical.Components.Substance O2(substanceData=Substances.Oxygen_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.000167)
          "free dissolved undound oxygen"
          annotation (Placement(transformation(extent={{96,20},{116,40}})));
        Chemical.Components.Membrane freeO2(KC=KC) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={118,0})));
        Chemical.Components.Substance O2_E(amountOfSubstance_start(displayUnit=
                "mmol") = 0.000125, substanceData=Substances.Oxygen_aqueous)
          "free dissolved undound O2"
          annotation (Placement(transformation(extent={{96,-38},{116,-18}})));
        Chemical.Components.Substance
                             K(substanceData=Substances.Potassium_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.004)
          annotation (Placement(transformation(extent={{-100,20},{-120,40}})));
        Chemical.Components.Substance
                             Na(                               substanceData=
              Substances.Sodium_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.138)
          annotation (Placement(transformation(extent={{-124,20},{-144,40}})));
        Chemical.Components.Substance
                             Na_E(substanceData=Substances.Sodium_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.007)
          annotation (Placement(transformation(extent={{-118,-38},{-138,-18}})));
        Chemical.Components.Substance
                             K_E(substanceData=Substances.Potassium_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.096)
          annotation (Placement(transformation(extent={{-112,-38},{-92,-18}})));
        Chemical.Components.Substance H2PO4_E(substanceData=Substances.DihydrogenPhosphate_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.000175)
          annotation (Placement(transformation(extent={{-84,-38},{-64,-18}})));
        Chemical.Components.Substance ADP_E(substanceData(z=-3),
            amountOfSubstance_start(displayUnit="mmol") = 9.6e-05)
          annotation (Placement(transformation(extent={{-114,-62},{-94,-42}})));
        Chemical.Components.Substance ATP_E(substanceData(
            z=-4,
            DfH_25degC=16700,
            DfG_25degC_1bar=30500,
            References={"http://www.wiley.com/college/pratt/0471393878/student/review/thermodynamics/7_relationship.html"}),
            amountOfSubstance_start(displayUnit="mmol") = 0.00128)
          annotation (Placement(transformation(extent={{-146,-62},{-166,-42}})));
        Chemical.Components.Substance HPO4_E(substanceData=Substances.HydrogenPhosphate_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.000495)
          annotation (Placement(transformation(extent={{-84,-62},{-64,-42}})));
        Chemical.Components.Substance H2PO4(substanceData=Substances.DihydrogenPhosphate_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.000365)
          annotation (Placement(transformation(extent={{-86,78},{-66,98}})));
        Chemical.Components.Substance HPO4(substanceData=Substances.HydrogenPhosphate_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.001635)
          annotation (Placement(transformation(extent={{-80,58},{-60,78}})));
        Chemical.Components.Substance globulins(substanceData(
            MolarWeight=34,
            z=-2.43,
            density=1080), amountOfSubstance_start(displayUnit="mmol") = 0.00082)
          annotation (Placement(transformation(extent={{150,76},{130,96}})));
        Chemical.Components.Substance Ca(substanceData=Substances.Calcium_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.00175) "Ca2+"
          annotation (Placement(transformation(extent={{-78,20},{-98,40}})));
        Chemical.Components.Substance Mg(substanceData=Substances.Magnesium_aqueous,
            amountOfSubstance_start(displayUnit="mmol") = 0.00108) "Mg2+"
          annotation (Placement(transformation(extent={{-112,-84},{-92,-64}})));
        Chemical.Components.Substance DPG(amountOfSubstance_start(displayUnit=
                "mmol") = 0.0051, substanceData(
            MolarWeight=0.266,
            z=-2.2,
            density=1000))
          annotation (Placement(transformation(extent={{128,-94},{108,-74}})));
        Chemical.Components.Substance GSH(substanceData(
            MolarWeight=0.2,
            z=-1,
            density=1000), amountOfSubstance_start(displayUnit="mmol") = 0.00223)
          annotation (Placement(transformation(extent={{164,-94},{144,-74}})));
        Components.Reaction          chemicalReaction1(nP=2) "10^(-6.66 + 3)"
          annotation (Placement(transformation(extent={{-56,78},{-36,98}})));
        Components.Reaction          HendersonHasselbalch(nP=2, nS=2,
        useKineticsInput=false) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-34,64},{-14,44}})));
        Sources.Buffer Hemoglobin(
          substanceData(z=1.045),
          a_start=10^(-7.2),
          BufferValue=3) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              origin={82,-76})));
      equation
      //  pH_p = -log10(H.a);
       // pH_e = -log10(H_E.a);
      connect(H2O.solution, blood_plasma.solution)
        annotation (Line(points={{-150,44},{-150,44},{-150,26},{-150,26},{-150,20},{108,
                20},{108,12.88}}, color={127,127,0}));
      connect(Cl.solution, blood_plasma.solution) annotation (Line(
          points={{0,20},{0,16},{0,12.88},{108,12.88}},
          color={127,127,0}));
        connect(H2O_E.solution, blood_erythrocytes.solution) annotation (Line(
              points={{-148,-38},{108,-38},{108,-99.1}},
                                                      color={127,127,0}));
        connect(Cl_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{0,-38},{0,-68},{0,-99.1},{108,-99.1}},
            color={127,127,0}));
        connect(HCO3_E.solution, blood_erythrocytes.solution) annotation (Line(
              points={{-24,-38},{108,-38},{108,-99.1}},
                                                    color={127,127,0}));
        connect(Aquapirin.port_b, H2O_E.port_a) annotation (Line(
            points={{-168,-10},{-168,-28},{-164,-28}},
            color={158,66,200},
            thickness=1));
        connect(Aquapirin.port_a, H2O.port_a) annotation (Line(
            points={{-168,10},{-168,54},{-166,54}},
            color={158,66,200},
            thickness=1));
        connect(Band3.port_a, HCO3.port_a) annotation (Line(
            points={{-6,10},{-6,30},{-8,30}},
            color={158,66,200},
            thickness=1));
        connect(Band3.port_b, HCO3_E.port_a) annotation (Line(
            points={{-6,-10},{-6,-28},{-8,-28}},
            color={158,66,200},
            thickness=1));
        connect(Band3_.port_b, Cl_E.port_a) annotation (Line(
            points={{18,-10},{18,-28},{16,-28}},
            color={158,66,200},
            thickness=1));
        connect(Band3_.port_a, Cl.port_a) annotation (Line(
            points={{18,10},{18,30},{16,30}},
            color={158,66,200},
            thickness=1));
      connect(HCO3.solution, blood_plasma.solution) annotation (Line(
          points={{-24,20},{108,20},{108,12.88}},
          color={127,127,0}));
        connect(blood_plasma.solution, permeableUncharged.solution) annotation (Line(
            points={{108,12.88},{108,20},{162,20}},
            color={127,127,0}));
        connect(blood_erythrocytes.solution, permeableUncharged_E.solution)
          annotation (Line(
            points={{108,-99.1},{108,-38},{160,-38}},
            color={127,127,0}));
        connect(blood_erythrocytes.solution,chargedImpermeable_E. solution)
          annotation (Line(
            points={{108,-99.1},{108,-38},{140,-38},{140,-62},{148,-62}},
            color={127,127,0}));
        connect(permeableUncharged.port_a, leak.port_a) annotation (Line(
            points={{146,30},{140,30},{140,10}},
            color={158,66,200},
            thickness=1));
        connect(permeableUncharged_E.port_a, leak.port_b) annotation (Line(
            points={{144,-28},{140,-28},{140,-10}},
            color={158,66,200},
            thickness=1));
        connect(MCT_.port_a, Lac.port_a) annotation (Line(
            points={{78,10},{78,30},{76,30}},
            color={158,66,200},
            thickness=1));
        connect(MCT_.port_b, Lac_E.port_a) annotation (Line(
            points={{78,-10},{78,-28},{76,-28}},
            color={158,66,200},
            thickness=1));
        connect(Lac.solution, blood_plasma.solution) annotation (Line(
            points={{60,20},{108,20},{108,12.88}},
            color={127,127,0}));
        connect(blood_erythrocytes.solution, Lac_E.solution) annotation (Line(
            points={{108,-99.1},{108,-38},{60,-38}},
            color={127,127,0}));
        connect(H_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{34,-38},{108,-38},{108,-99.1}},
            color={127,127,0}));
        connect(H_E.port_a, MCT.port_b) annotation (Line(
            points={{50,-28},{52,-28},{52,-10}},
            color={158,66,200},
            thickness=1));
        connect(MCT.port_a, H.port_a) annotation (Line(
            points={{52,10},{52,30},{50,30}},
            color={158,66,200},
            thickness=1));
        connect(blood_plasma.solution, H.solution) annotation (Line(
            points={{108,12.88},{108,20},{34,20}},
            color={127,127,0}));
        connect(CO2.port_a, freeCO2.port_a) annotation (Line(
            points={{-40,30},{-38,30},{-38,12}},
            color={158,66,200},
            thickness=1));
        connect(freeCO2.port_b, CO2_E.port_a) annotation (Line(
            points={{-38,-8},{-38,-28}},
            color={158,66,200},
            thickness=1));
        connect(blood_plasma.solution, CO2.solution) annotation (Line(
            points={{108,12.88},{108,20},{-56,20}},
            color={127,127,0}));
        connect(CO2_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{-54,-38},{108,-38},{108,-99.1}},
            color={127,127,0}));
        connect(blood_plasma.solution, O2.solution) annotation (Line(
            points={{108,12.88},{108,20},{100,20}},
            color={127,127,0}));
        connect(O2_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{100,-38},{108,-38},{108,-99.1}},
            color={127,127,0}));
        connect(O2_E.port_a, freeO2.port_b) annotation (Line(
            points={{116,-28},{118,-28},{118,-10}},
            color={158,66,200},
            thickness=1));
        connect(freeO2.port_a, O2.port_a) annotation (Line(
            points={{118,10},{118,30},{116,30}},
            color={158,66,200},
            thickness=1));
        connect(H2O.solution, K.solution) annotation (Line(
            points={{-150,44},{-150,20},{-104,20}},
            color={158,66,200}));
        connect(H2O.solution, Na.solution) annotation (Line(
            points={{-150,44},{-150,20},{-128,20}},
            color={158,66,200}));
        connect(H2O_E.solution, Na_E.solution) annotation (Line(
            points={{-148,-38},{-122,-38}},
            color={158,66,200}));
        connect(H2O_E.solution, K_E.solution) annotation (Line(
            points={{-148,-38},{-108,-38}},
            color={158,66,200}));
        connect(H2O_E.solution, H2PO4_E.solution) annotation (Line(
            points={{-148,-38},{-80,-38}},
            color={127,127,0}));
        connect(ADP_E.solution, K_E.solution) annotation (Line(
            points={{-110,-62},{-110,-38},{-108,-38}},
            color={158,66,200}));
        connect(ATP_E.solution, Na_E.solution) annotation (Line(
            points={{-150,-62},{-122,-62},{-122,-38}},
            color={127,127,0}));
        connect(H2O.solution, H2PO4.solution) annotation (Line(
            points={{-150,44},{-150,44},{-150,20},{-82,20},{-82,78}},
            color={127,127,0}));
        connect(HPO4_E.solution, H2PO4_E.solution) annotation (Line(
            points={{-80,-62},{-110,-62},{-110,-38},{-80,-38}},
            color={127,127,0}));
        connect(HPO4.solution, H2PO4.solution) annotation (Line(
            points={{-76,58},{-82,58},{-82,78}},
            color={127,127,0}));
        connect(Ca.solution, CO2.solution) annotation (Line(
            points={{-82,20},{-82,20},{-56,20}},
            color={127,127,0}));
        connect(Mg.solution, blood_erythrocytes.solution) annotation (Line(
            points={{-108,-84},{-108,-38},{108,-38},{108,-99.1}},
            color={127,127,0}));
        connect(DPG.solution, permeableUncharged_E.solution) annotation (Line(
            points={{124,-94},{140,-94},{140,-38},{160,-38}},
            color={127,127,0}));
        connect(GSH.solution, permeableUncharged_E.solution) annotation (Line(
            points={{160,-94},{140,-94},{140,-38},{160,-38}},
            color={127,127,0}));
        connect(H2PO4.port_a, chemicalReaction1.substrates[1]) annotation (Line(
            points={{-66,88},{-58,88},{-56,88}},
            color={158,66,200},
            thickness=1));
        connect(HPO4.port_a, chemicalReaction1.products[1]) annotation (Line(
            points={{-60,68},{-60,68},{-22,68},{-22,86},{-36,86}},
            color={158,66,200},
            thickness=1));
        connect(H.port_a, chemicalReaction1.products[2]) annotation (Line(
            points={{50,30},{50,30},{50,90},{50,90},{-36,90}},
            color={158,66,200},
            thickness=1));
        connect(CO2.port_a, HendersonHasselbalch.substrates[2]) annotation (Line(
            points={{-40,30},{-38,30},{-38,52},{-34,52}},
            color={158,66,200},
            thickness=1));
        connect(HCO3.port_a, HendersonHasselbalch.products[2]) annotation (Line(
            points={{-8,30},{-8,30},{-8,52},{-14,52}},
            color={158,66,200},
            thickness=1));
        connect(HendersonHasselbalch.substrates[1], H2O.port_a) annotation (Line(
            points={{-34,56},{-166,56},{-166,54}},
            color={158,66,200},
            thickness=1));
        connect(HendersonHasselbalch.products[1], H.port_a) annotation (Line(
            points={{-14,56},{18,56},{50,56},{50,30}},
            color={158,66,200},
            thickness=1));
        connect(Hemoglobin.solution, blood_erythrocytes.solution) annotation (Line(
              points={{88,-86},{94,-86},{108,-86},{108,-99.1}}, color={127,127,
                0}));
        connect(Hemoglobin.port_a, H_E.port_a) annotation (Line(points={{72,-76},{64,-76},
                {50,-76},{50,-28}}, color={158,66,200}));

        connect(albumin.solution, blood_plasma.solution) annotation (Line(
            points={{108,76},{126,76},{126,20},{108,20},{108,12.88}},
            color={127,127,0},
            smooth=Smooth.None));
        connect(globulins.solution, blood_plasma.solution) annotation (Line(
            points={{146,76},{126,76},{126,20},{108,20},{108,12.88}},
            color={127,127,0},
            smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>Blood eqiulibrium across erythrocyte membrane bewteen blood plasma and intracellular fluid of erythrocytes.</p>
<p>Data of blood status are from:</p>
<p>Raftos, J.E., Bulliman, B.T. and Kuchel, P.W. Evaluation of an electrochemical model of erythrocyte pH buffering using 31P nuclear magnetic resonance data. <i>The Journal of general physiology</i> 1990;95(6):1183-1204. </p>
</html>",      revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=1e-008),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},
                  {180,100}}), graphics));
      end RedCellMembrane;
    end Dev;
  end AcidBase;

  package Hemoglobin "Hemoglobin blood gases binding"
    model Allosteric_Hemoglobin_MWC "Monod,Wyman,Changeux (1965)"
      extends Modelica.Icons.Example;

      constant Modelica.SIunits.AmountOfSubstance THb = 0.001
        "Total amount of hemoglobin";

      constant Modelica.SIunits.Temperature T=298.15 "Base Temperature";
      constant Real RT=Modelica.Constants.R*T;

      constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        "Amount of solution used for molarity to mole fraction conversion";

      constant Modelica.SIunits.Volume OneLiter = 0.001;

      constant Real L=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      constant Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      constant Modelica.SIunits.Concentration KR=0.000671946
        "Oxygen dissociation coefficient on relaxed(R) hemoglobin subunit";

      constant Real KRx=KR*OneLiter/AmountOfSolutionIn1L
        "Mole fraction based KR";

    //Relative Gibbs formation energies of the substances in the system:
      constant Modelica.SIunits.MolarEnergy
        GO2aq=-RT*log(0.0013/55.508),
        GR0=0,                            GT0=GR0 -RT*log(L),
        GR1=GR0+GO2aq +RT*log(KRx/4),     GT1=GR1 -RT*log(c*L),
        GR2=GR1+GO2aq +RT*log(KRx/(3/2)), GT2=GR2 -RT*log(c^2*L),
        GR3=GR2+GO2aq +RT*log(KRx/(2/3)), GT3=GR3 -RT*log(c^3*L),
        GR4=GR3+GO2aq +RT*log(KRx*4),     GT4=GR4 -RT*log(c^4*L);
                                      //*0.018),

      parameter Real KC = 0.001 "Slow down factor";

      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-72,-102},{94,124}})));

      Chemical.Components.Substance oxygen_unbound(substanceData(
            DfG_25degC_1bar=GO2aq), amountOfSubstance_start(displayUnit="mol")=
             1e-5)
        annotation (Placement(transformation(extent={{-62,-46},{-42,-26}})));

      Chemical.Components.Substance T0(substanceData(DfG_25degC_1bar=GT0),
          amountOfSubstance_start=THb)
        annotation (Placement(transformation(extent={{34,78},{54,98}})));

      Chemical.Components.Substance T1(substanceData(DfG_25degC_1bar=GT1),
          amountOfSubstance_start=THb*1e-4)
        annotation (Placement(transformation(extent={{34,36},{54,56}})));

      Chemical.Components.Substance T2(substanceData(DfG_25degC_1bar=GT2),
          amountOfSubstance_start=THb*1e-8)
        annotation (Placement(transformation(extent={{34,-10},{54,10}})));

      Chemical.Components.Substance R1(substanceData(DfG_25degC_1bar=GR1),
          amountOfSubstance_start=THb*1e-8)
        annotation (Placement(transformation(extent={{-20,36},{0,56}})));

      Chemical.Components.Substance R2(substanceData(DfG_25degC_1bar=GR2),
          amountOfSubstance_start=THb*1e-10)
        annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

      Chemical.Components.Substance T3(substanceData(DfG_25degC_1bar=GT3),
          amountOfSubstance_start=THb*1e-12)
        annotation (Placement(transformation(extent={{34,-54},{54,-34}})));

      Chemical.Components.Substance R3(substanceData(DfG_25degC_1bar=GR3),
          amountOfSubstance_start=THb*1e-12)
        annotation (Placement(transformation(extent={{-20,-54},{0,-34}})));

      Chemical.Components.Substance T4(substanceData(DfG_25degC_1bar=GT4),
          amountOfSubstance_start=THb*1e-17)
        annotation (Placement(transformation(extent={{34,-92},{54,-72}})));

      Chemical.Components.Substance R4(substanceData(DfG_25degC_1bar=GR4),
          amountOfSubstance_start=THb*1e-14)
        annotation (Placement(transformation(extent={{-20,-92},{0,-72}})));

      Chemical.Components.Substance R0(substanceData(DfG_25degC_1bar=GR0),
          amountOfSubstance_start=THb*1e-7)
        annotation (Placement(transformation(extent={{-20,78},{0,98}})));

      Chemical.Components.Reaction quaternaryForm(KC=KC)
        annotation (Placement(transformation(extent={{4,78},{24,98}})));
      Chemical.Components.Reaction oxyR1(nP=2, KC=KC) annotation (Placement(
            transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-8,64})));
      Chemical.Components.Reaction oxyT1(nP=2, KC=KC) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={44,64})));
      Chemical.Components.Reaction oxyR2(nP=2, KC=KC) annotation (Placement(
            transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-10,22})));
      Chemical.Components.Reaction oxyR3(nP=2, KC=KC) annotation (Placement(
            transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-10,-24})));
      Chemical.Components.Reaction oxyR4(nP=2, KC=KC) annotation (Placement(
            transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-10,-66})));
      Chemical.Components.Reaction oxyT2(nP=2, KC=KC) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={44,22})));
      Chemical.Components.Reaction oxyT3(nP=2, KC=KC) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={44,-24})));
      Chemical.Components.Reaction oxyT4(nP=2, KC=KC) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={44,-66})));
      Chemical.Components.Reaction quaternaryForm1(KC=KC)
        annotation (Placement(transformation(extent={{8,36},{28,56}})));
      Chemical.Components.Reaction quaternaryForm2(KC=KC)
        annotation (Placement(transformation(extent={{8,-10},{28,10}})));
      Chemical.Components.Reaction quaternaryForm3(KC=KC)
        annotation (Placement(transformation(extent={{8,-54},{28,-34}})));
      Chemical.Components.Reaction quaternaryForm4(KC=KC)
        annotation (Placement(transformation(extent={{10,-92},{30,-72}})));

      Modelica.Blocks.Sources.Clock clock(offset=10)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-84,62})));
    Chemical.Sources.ExternalIdealGasSubstance O2_in_air(
        TotalPressure(displayUnit="kPa") = 101325.0144354,
        substanceData=Substances.Oxygen_gas,
        PartialPressure(displayUnit="kPa") = 1000,
        usePartialPressureInput=true) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-84,22})));

      Chemical.Components.GasSolubility gasSolubility(useWaterCorrection=false,
          KC=KC)
        annotation (Placement(transformation(extent={{-94,-16},{-74,4}})));

      Chemical.Components.Substance H2O(substanceData=Substances.Water_liquid,
          amountOfSubstance_start=38.7)
        annotation (Placement(transformation(extent={{64,-100},{84,-80}})));

      Real sO2;
    equation
      sO2 = (R1.x + 2*R2.x + 3*R3.x + 4*R4.x + T1.x + 2*T2.x + 3*T3.x + 4*T4.x) /
       (4*(R0.x + R1.x + R2.x + R3.x + R4.x + T0.x + T1.x + T2.x + T3.x + T4.x));

      connect(quaternaryForm.products[1],T0. port_a) annotation (Line(
          points={{24,88},{54,88}},
          color={107,45,134}));
      connect(oxyR1.substrates[1],R1. port_a) annotation (Line(
          points={{-8,54},{-8,46},{0,46}},
          color={107,45,134}));
      connect(R1.port_a,oxyR2. products[1]) annotation (Line(
          points={{0,46},{0,32},{-12,32}},
          color={107,45,134}));
      connect(oxyR2.substrates[1],R2. port_a) annotation (Line(
          points={{-10,12},{-10,0},{0,0}},
          color={107,45,134}));
      connect(oxyR3.substrates[1],R3. port_a) annotation (Line(
          points={{-10,-34},{-10,-44},{0,-44}},
          color={107,45,134}));
      connect(oxyR3.products[1],R2. port_a) annotation (Line(
          points={{-12,-14},{-12,-7},{0,-7},{0,0}},
          color={107,45,134}));
      connect(R3.port_a,oxyR4. products[1]) annotation (Line(
          points={{0,-44},{0,-56},{-12,-56}},
          color={107,45,134}));
      connect(oxyR4.substrates[1],R4. port_a) annotation (Line(
          points={{-10,-76},{-10,-82},{0,-82}},
          color={107,45,134}));
      connect(oxyT1.products[1],T0. port_a) annotation (Line(
          points={{46,74},{46,88},{54,88}},
          color={107,45,134}));
      connect(oxyT1.substrates[1],T1. port_a) annotation (Line(
          points={{44,54},{44,46},{54,46}},
          color={107,45,134}));
      connect(T1.port_a,oxyT2. products[1]) annotation (Line(
          points={{54,46},{54,32},{46,32}},
          color={107,45,134}));
      connect(oxyT3.substrates[1],T3. port_a) annotation (Line(
          points={{44,-34},{44,-44},{54,-44}},
          color={107,45,134}));
      connect(T3.port_a,oxyT4. products[1]) annotation (Line(
          points={{54,-44},{54,-56},{46,-56}},
          color={107,45,134}));
      connect(oxyT4.substrates[1],T4. port_a) annotation (Line(
          points={{44,-76},{44,-82},{54,-82}},
          color={107,45,134}));
      connect(R0.port_a,quaternaryForm. substrates[1]) annotation (Line(
          points={{0,88},{4,88}},
          color={107,45,134}));
      connect(R0.port_a,oxyR1. products[1]) annotation (Line(
          points={{0,88},{0,74},{-10,74}},
          color={107,45,134}));
      connect(R1.port_a,quaternaryForm1. substrates[1]) annotation (Line(
          points={{0,46},{8,46}},
          color={107,45,134}));
      connect(quaternaryForm1.products[1],T1. port_a) annotation (Line(
          points={{28,46},{54,46}},
          color={107,45,134}));
      connect(R2.port_a,quaternaryForm2. substrates[1]) annotation (Line(
          points={{0,0},{8,0}},
          color={107,45,134}));
      connect(R3.port_a,quaternaryForm3. substrates[1]) annotation (Line(
          points={{0,-44},{8,-44}},
          color={107,45,134}));
      connect(quaternaryForm3.products[1],T3. port_a) annotation (Line(
          points={{28,-44},{54,-44}},
          color={107,45,134}));
      connect(R4.port_a,quaternaryForm4. substrates[1]) annotation (Line(
          points={{0,-82},{10,-82}},
          color={107,45,134}));
      connect(quaternaryForm4.products[1],T4. port_a) annotation (Line(
          points={{30,-82},{54,-82}},
          color={107,45,134}));
      connect(oxyR1.products[2],oxygen_unbound. port_a)
                                          annotation (Line(
          points={{-6,74},{-42,74},{-42,-36}},
          color={107,45,134}));
      connect(oxyR2.products[2],oxygen_unbound. port_a)
                                          annotation (Line(
          points={{-8,32},{-42,32},{-42,-36}},
          color={107,45,134}));
      connect(oxyR3.products[2],oxygen_unbound. port_a)
                                          annotation (Line(
          points={{-8,-14},{-42,-14},{-42,-36}},
          color={107,45,134}));
      connect(oxyR4.products[2],oxygen_unbound. port_a)
                                          annotation (Line(
          points={{-8,-56},{-42,-56},{-42,-36}},
          color={107,45,134}));
      connect(oxygen_unbound.port_a, oxyT1.products[2])
                                          annotation (Line(
          points={{-42,-36},{-42,74},{42,74}},
          color={107,45,134}));
      connect(oxygen_unbound.port_a, oxyT2.products[2])
                                          annotation (Line(
          points={{-42,-36},{-42,32},{42,32}},
          color={107,45,134}));
      connect(oxygen_unbound.port_a, oxyT3.products[2])
                                          annotation (Line(
          points={{-42,-36},{-42,-14},{42,-14}},
          color={107,45,134}));
      connect(oxygen_unbound.port_a, oxyT4.products[2])
                                          annotation (Line(
          points={{-42,-36},{-42,-56},{42,-56}},
          color={107,45,134}));
      connect(O2_in_air.port_a, gasSolubility.gas_port) annotation (Line(
          points={{-84,12},{-84,4}},
          color={158,66,200}));
      connect(gasSolubility.liquid_port, oxygen_unbound.port_a) annotation (Line(
          points={{-84,-16},{-84,-36},{-42,-36}},
          color={158,66,200}));
      connect(oxygen_unbound.solution, solution.solution) annotation (Line(
          points={{-58,-46},{-58,-64},{-58,-64},{-58,-82},{-16,-82},{-16,-100},
              {-16,-100},{-16,-100},{-16,-99.74},{22,-99.74},{60.8,-99.74}},
          color={127,127,0}));
      connect(R0.solution, solution.solution) annotation (Line(
          points={{-16,78},{-16,-99.74},{60.8,-99.74}},
          color={127,127,0}));
      connect(T0.solution, solution.solution) annotation (Line(
          points={{38,78},{38,-99.74},{60.8,-99.74}},
          color={127,127,0}));
      connect(R1.solution, solution.solution) annotation (Line(points={{-16,36},
              {-16,-99.74},{60.8,-99.74}},
                                color={127,127,0}));
      connect(T1.solution, solution.solution) annotation (Line(points={{38,36},
              {38,-99.74},{60.8,-99.74}},
                          color={127,127,0}));
      connect(R2.solution, solution.solution) annotation (Line(points={{-16,-10},
              {-16,-99.74},{60.8,-99.74}},
                                color={127,127,0}));
      connect(T3.solution, solution.solution) annotation (Line(points={{38,-54},
              {38,-99.74},{60.8,-99.74}},
                                color={127,127,0}));
      connect(R3.solution, solution.solution) annotation (Line(points={{-16,-54},
              {-16,-99.74},{60.8,-99.74}},
                                color={127,127,0}));
      connect(R4.solution, solution.solution) annotation (Line(points={{-16,-92},
              {-16,-99.74},{60.8,-99.74}},
                                color={127,127,0}));
      connect(T4.solution, solution.solution) annotation (Line(points={{38,-92},
              {38,-99.74},{60.8,-99.74}},
                                color={127,127,0}));
      connect(quaternaryForm2.products[1], T2.port_a) annotation (Line(
          points={{28,0},{54,0}},
          color={158,66,200}));
      connect(oxyT2.substrates[1], T2.port_a) annotation (Line(
          points={{44,12},{44,0},{54,0}},
          color={158,66,200}));
      connect(T2.port_a, oxyT3.products[1]) annotation (Line(
          points={{54,0},{54,-14},{46,-14}},
          color={158,66,200}));
      connect(T2.solution, solution.solution) annotation (Line(points={{38,-10},
              {38,-99.74},{60.8,-99.74}},
                                color={127,127,0}));
      connect(clock.y, O2_in_air.partialPressure) annotation (Line(
          points={{-84,51},{-84,32}},
          color={0,0,127}));
    connect(H2O.solution, solution.solution) annotation (Line(
        points={{68,-100},{68,-99.74},{60.8,-99.74}},
        color={127,127,0}));
      annotation (          experiment(StopTime=15000, Tolerance=0.01),
                              Documentation(info="<html>
<p>To understand the model is necessary to study the principles of MWC allosteric transitions first published by </p>
<p>[1] Monod,Wyman,Changeux (1965). &quot;On the nature of allosteric transitions: a plausible model.&quot; Journal of molecular biology 12(1): 88-118.</p>
<p><br>In short it is about binding oxygen to hemoglobin.</p>
<p>Oxgen are driven by its partial pressure using clock source - from very little pressure to pressure of 10kPa.</p>
<p>(Partial pressure of oxygen in air is the air pressure multiplied by the fraction of the oxygen in air.)</p>
<p>Hemoglobin was observed (by Perutz) in two structuraly different forms R and T.</p>
<p>These forms are represented by blocks T0..T4 and R0..R4, where the suffexed index means the number of oxygen bounded to the form.</p>
<p><br>In equilibrated model can be four chemical reactions removed and the results will be the same, but dynamics will change a lot. ;)</p>
<p>If you remove the quaternaryForm1,quaternaryForm2,quaternaryForm3,quaternaryForm4 then the model in equilibrium will be exactly the same as in MWC article.</p>
<p><br>Parameters was fitted to data of Severinghaus article from 1979. (For example at pO2=26mmHg is oxygen saturation sO2 = 48.27 %).</p>
</html>", revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Allosteric_Hemoglobin_MWC;

    model Allosteric_Hemoglobin2_MWC
      "Monod,Wyman,Changeux (1965) - The same allosteric hemoglobin model as Allosteric_Hemoglobin_MWC implemented by Speciation blocks"
      extends Modelica.Icons.Example;

      constant Modelica.SIunits.AmountOfSubstance THb = 0.001
        "Total amount of hemoglobin";

      constant Modelica.SIunits.Temperature T=298.15 "Base Temperature";
      constant Real RT=Modelica.Constants.R*T;

      constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        "Amount of solution used for molarity to mole fraction conversion";

      constant Modelica.SIunits.Volume OneLiter = 0.001;

      parameter Real L=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.SIunits.Concentration KR=0.000671946
        "oxygen dissociation on relaxed(R) hemoglobin subunit";
                                                                  //*7.875647668393782383419689119171e-5
                                                                //10.500001495896 7.8756465463794e-05

      parameter Modelica.SIunits.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.SIunits.MoleFraction KRx = KR*OneLiter/AmountOfSolutionIn1L;
      parameter Modelica.SIunits.MoleFraction KTx = KT*OneLiter/AmountOfSolutionIn1L;

      parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508);
      parameter Modelica.SIunits.ChemicalPotential DfG_uR = 0;
      parameter Modelica.SIunits.ChemicalPotential DfG_uRO2 = DfG_uR + DfG_O2 + RT * log(KRx);
      parameter Modelica.SIunits.ChemicalPotential DfG_uT = 0;
      parameter Modelica.SIunits.ChemicalPotential DfG_uTO2 = DfG_uT + DfG_O2 + RT * log(KTx);
      parameter Modelica.SIunits.ChemicalPotential DfG_tT = 0;
      parameter Modelica.SIunits.ChemicalPotential DfG_tR = DfG_tT + RT * log(L);

      parameter Real KC = 1e-3 "Slow down factor";
                               //0.000001

      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,42}})));

      Chemical.Components.Reaction quaternaryForm(KC=KC)
        annotation (Placement(transformation(extent={{12,-58},{32,-38}})));
      Chemical.Components.Speciation R0_in_R(NumberOfSubunits=4)
        annotation (Placement(transformation(extent={{-46,-48},{-26,-28}})));
       // AmountOfSubstance_start=4e-11)
      Chemical.Components.Speciation T0_in_T(NumberOfSubunits=4)
        annotation (Placement(transformation(extent={{76,-48},{56,-28}})));
       // AmountOfSubstance_start=totalAmountOfHemoglobin)
      Chemical.Components.Substance OxyRHm[4](
        each amountOfSubstance_start=5.88e-9,
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_O2 + RT*log(KRx) + DfG_tR/4))
        "Oxygenated subunit in R structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{-96,-18},{-76,2}})));

      Chemical.Components.Reaction oxygenation_R[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{-68,-18},{-48,2}})));
      Chemical.Components.Substance DeoxyRHm[4](
        each amountOfSubstance_start=1.58e-7,
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_tR/4))
        "Deoxygenated subunit in R structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{-8,-20},{-28,0}})));

      Chemical.Components.Substance OxyTHm[4](
        each amountOfSubstance_start=1e-4,
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_O2 + RT*log(KTx) + DfG_tT/4))
        "Oxygenated subunit in T structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{14,-18},{34,2}})));

      Chemical.Components.Reaction oxygenation_T[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{42,-18},{62,2}})));
      Chemical.Components.Substance DeoxyTHm[4](
        each amountOfSubstance_start=THb - 1e-4 - 1.58e-7 - 5.88e-9,
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_tT/4))
        "Deoxygenated subunit in T structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{96,-20},{76,0}})));

      Chemical.Components.Substance oxygen_unbound(substanceData(
            DfG_25degC_1bar=DfG_O2), amountOfSubstance_start=2e-8)
        annotation (Placement(transformation(extent={{-2,6},{18,26}})));
      Modelica.Blocks.Sources.Clock clock(offset=100)
        annotation (Placement(transformation(extent={{-40,74},{-20,94}})));
    Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
          usePartialPressureInput=true, substanceData=Substances.Oxygen_gas)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={8,68})));
      Chemical.Components.GasSolubility partialPressure1(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={8,40})));

      Real sO2 "Hemoglobin oxygen saturation";
      Chemical.Components.Substance H2O(substanceData=Substances.Water_liquid,
          amountOfSubstance_start=38.7)
        annotation (Placement(transformation(extent={{32,-92},{52,-72}})));
    equation
      sO2 = (sum(OxyRHm.x) + sum(OxyTHm.x)) /
      (sum(DeoxyRHm.x) + sum(DeoxyTHm.x) + sum(OxyRHm.x) + sum(OxyTHm.x));

      connect(OxyTHm.port_a, oxygenation_T.substrates[1])
                                               annotation (Line(
          points={{34,-8},{42,-8}},
          color={107,45,134}));
      connect(oxygenation_T.products[1], DeoxyTHm.port_a)
                                             annotation (Line(
          points={{62,-10},{70,-10},{76,-10}},
          color={107,45,134}));

      connect(clock.y, oxygen_in_air.partialPressure) annotation (Line(
          points={{-19,84},{8,84},{8,78}},
          color={0,0,127}));
      connect(OxyRHm.port_a, oxygenation_R.substrates[1]) annotation (Line(
          points={{-76,-8},{-68,-8}},
          color={107,45,134}));
      connect(DeoxyRHm.port_a, R0_in_R.subunits) annotation (Line(
          points={{-28,-10},{-39,-10},{-39,-27.8}},
          color={107,45,134}));
      connect(oxygenation_R.products[1], DeoxyRHm.port_a) annotation (Line(
          points={{-48,-10},{-28,-10}},
          color={107,45,134}));

      connect(T0_in_T.subunits, DeoxyTHm.port_a)   annotation (Line(
          points={{69,-27.8},{69,-10},{76,-10}},
          color={107,45,134}));

      connect(oxygen_in_air.port_a, partialPressure1.gas_port) annotation (Line(
          points={{8,58},{8,50}},
          color={158,66,200}));
      connect(partialPressure1.liquid_port, oxygen_unbound.port_a) annotation (Line(
          points={{8,30},{8,16},{18,16}},
          color={158,66,200}));
      connect(R0_in_R.port_a, quaternaryForm.substrates[1]) annotation (Line(
          points={{-26,-48},{-26,-48},{12,-48}},
          color={158,66,200}));
      connect(quaternaryForm.products[1], T0_in_T.port_a) annotation (Line(
          points={{32,-48},{32,-48},{56,-48}},
          color={158,66,200}));

      for i in 1:4 loop
        connect(oxygenation_T[i].products[2], oxygen_unbound.port_a) annotation (Line(
          points={{62,-6},{70,-6},{70,16},{18,16}},
          color={107,45,134}));
        connect(oxygenation_R[i].products[2], oxygen_unbound.port_a) annotation (Line(
          points={{-48,-6},{-34,-6},{-34,16},{18,16}},
          color={107,45,134}));
      connect(R0_in_R.subunitSolution, DeoxyRHm[i].solution) annotation (Line(
          points={{-32,-32},{-32,-22},{-12,-22},{-12,-20}},
          color={127,127,0}));
      connect(R0_in_R.subunitSolution, OxyRHm[i].solution) annotation (Line(
          points={{-32,-32},{-32,-22},{-92,-22},{-92,-18}},
          color={127,127,0}));
      connect(OxyTHm[i].solution, T0_in_T.subunitSolution) annotation (Line(
          points={{18,-18},{18,-22},{62,-22},{62,-32}},
          color={127,127,0}));
      connect(DeoxyTHm[i].solution, T0_in_T.subunitSolution) annotation (Line(
          points={{92,-20},{92,-22},{62,-22},{62,-32}},
          color={127,127,0}));
      end for;

      connect(R0_in_R.solution, solution.solution) annotation (Line(
          points={{-42,-48},{-42,-98.58},{60,-98.58}},
          color={127,127,0}));
      connect(T0_in_T.solution, solution.solution) annotation (Line(
          points={{72,-48},{72,-98.58},{60,-98.58}},
          color={127,127,0}));
      connect(oxygen_unbound.solution, solution.solution) annotation (Line(points={{2,6},{2,
              -98.58},{60,-98.58}},            color={127,127,0}));
      connect(solution.solution, H2O.solution) annotation (Line(
          points={{60,-98.58},{36,-98.58},{36,-92}},
          color={127,127,0}));

      annotation (          experiment(StopTime=15000),
        Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<p><br>To understand the model is necessary to study the principles of MWC allosteric transitions first published by </p>
<p>[1] Monod,Wyman,Changeux (1965). &quot;On the nature of allosteric transitions: a plausible model.&quot; Journal of molecular biology 12(1): 88-118.</p>
<p><br>In short it is about binding oxygen to hemoglobin.</p>
<p>Oxgen are driven by its partial pressure using clock source - from very little pressure to pressure of 10kPa.</p>
<p>(Partial pressure of oxygen in air is the air pressure multiplied by the fraction of the oxygen in air.)</p>
<p>Hemoglobin was observed (by Perutz) in two structuraly different forms R and T.</p>
<p>These forms are represented by blocks T0..T4 and R0..R4, where the suffexed index means the number of oxygen bounded to the form.</p>
<p><br>In equilibrated model can be four chemical reactions removed and the results will be the same, but dynamics will change a lot. ;)</p>
<p>If you remove the quaternaryForm1,quaternaryForm2,quaternaryForm3,quaternaryForm4 then the model in equilibrium will be exactly the same as in MWC article.</p>
<p><br>Parameters was fitted to data of Severinghaus article from 1979. (For example at pO2=26mmHg is oxygen saturation sO2 = 48.27 %).</p>
</html>"));
    end Allosteric_Hemoglobin2_MWC;

    model HemoglobinQuaternaryForm
      "Hemoglobib quaternary form - part of multiple-ligand allosteric hemoglobin model"

      constant Integer N=12
        "Number of distinguished independent sides in quaternary structure";
      constant Real RT=Modelica.Constants.R*298.15;

      parameter Modelica.SIunits.MolarEnthalpy Ho = 59000
        "Enthalpy of deoxygenation";
      parameter Modelica.SIunits.MoleFraction Ko37 "KRx and KTx at 37degC";
      parameter Modelica.SIunits.MoleFraction Ko25 = Ko37*exp((Ho/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRx and KTx at 25degC";

      parameter Modelica.SIunits.MolarEnthalpy Hco = 59000
        "Enthalpy of carbon monoxide dissociation";
      parameter Modelica.SIunits.MoleFraction Kco37
        "Carboxyhemoglobin dissociation at 37degC";
      parameter Modelica.SIunits.MoleFraction Kco25 = Kco37*exp((Hco/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "Carboxyhemoglobin dissociation at 25degC";

      parameter Modelica.SIunits.MolarEnthalpy Hh
        "Enthalpy of deprotonation of h site";
      parameter Modelica.SIunits.MoleFraction Kh37 "KRhx and KThx at 37 degC";
      parameter Modelica.SIunits.MoleFraction Kh25 = Kh37*exp(((Hh)/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRhx and KThx at 25 degC";

      parameter Modelica.SIunits.MolarEnthalpy Hz
        "Enthalpy of deprotonation of -NH3+ terminus";
      parameter Modelica.SIunits.MoleFraction Kz37 "KRzx and KTzx at 37 degC";
      parameter Modelica.SIunits.MoleFraction Kz25 = Kz37*exp(((Hz)/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRzx and KTzx at 25 degC";

      parameter Modelica.SIunits.MolarEnthalpy Hc "Enthalpy of carboxylation";
      parameter Modelica.SIunits.MoleFraction Kc37 "KRcx and KTcx at 37degC";
      parameter Modelica.SIunits.MoleFraction Kc25 = Kc37*exp((Hc/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRcx and KTcx at 25degC";

      parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508)  + Chemical.Examples.Substances.Oxygen_gas.DfG_25degC_1bar;
      parameter Modelica.SIunits.ChemicalPotential DfH_O2 = Chemical.Examples.Substances.Oxygen_aqueous.DfH_25degC;

      parameter Modelica.SIunits.ChemicalPotential DfG_CO2 = -RT*log(0.034/55.508)  + Chemical.Examples.Substances.CarbonDioxide_gas.DfG_25degC_1bar;
      parameter Modelica.SIunits.ChemicalPotential DfH_CO2 = Chemical.Examples.Substances.CarbonDioxide_aqueous.DfH_25degC;

      parameter Modelica.SIunits.ChemicalPotential DfG_selectedForm
        "DfG_tR and DfG_tT";
      parameter Modelica.SIunits.MolarEnthalpy DfH_selectedForm = 0
        "DfH_tR and DfH_tT";

      parameter Real KC = 1e-3 "Slow down factor";
                               //0.000001
      parameter Modelica.SIunits.MoleFraction initialO2
        "Initial mole fraction of unbound oxygen disoluted around hemoglobin";
      parameter Modelica.SIunits.MoleFraction initialH
        "Initial mole fraction of H+";
      parameter Modelica.SIunits.MoleFraction initialCO2
        "Initial mole fraction of unbound carbon dioxide disoluted around hemoglobin";
      parameter Modelica.SIunits.AmountOfSubstance initialHb
        "Initial amount of hemoglobin tetramers in this quaternary form";

      Chemical.Components.Speciation speciation(NumberOfSubunits=N)
        annotation (Placement(transformation(extent={{-18,-72},{2,-52}})));
       // AmountOfSubstance_start=4e-11)
       // AmountOfSubstance_start=totalAmountOfHemoglobin)
      Chemical.Components.Substance OxyHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_O2 + RT*log(Ko25) + DfG_selectedForm/N,
              DfH_25degC=DfH_O2 -Ho + DfH_selectedForm/N),
        each amountOfSubstance_start=initialO2*initialHb/(Ko37 +
            initialO2)) "Oxygenated subunit"
        annotation (Placement(transformation(extent={{-88,14},{-68,34}})));

      Chemical.Components.Reaction o[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
      Chemical.Components.Substance DeoxyHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_selectedForm/N,
              DfH_25degC=DfH_selectedForm/N),
        each amountOfSubstance_start=Ko37*initialHb/(Ko37 + initialO2))
        "Deoxygenated subunit"
        annotation (Placement(transformation(extent={{-8,12},{-28,32}})));

      Chemical.Interfaces.SolutionPort solution
        annotation (Placement(transformation(extent={{-50,-82},{-32,-62}}),
            iconTransformation(extent={{-50,-90},{-30,-70}})));
      Chemical.Interfaces.SubstancePort_b O2 annotation (Placement(transformation(
              extent={{-28,32},{-8,52}}),  iconTransformation(extent={{-90,70},{-70,
                90}})));
      Chemical.Interfaces.SubstancePort_a selectedForm
        annotation (Placement(transformation(extent={{26,-82},{46,-62}}),
            iconTransformation(extent={{30,-90},{50,-70}})));
      Chemical.Components.Substance HmAH[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=RT*log(Kh25) + DfG_selectedForm/N,
              DfH_25degC=-Hh+DfH_selectedForm/N),
        each amountOfSubstance_start=initialH*initialHb/(Kh37 + initialH))
        "Protonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{10,12},{30,32}})));
      Chemical.Components.Reaction h[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{36,32},{56,12}})));
      Chemical.Components.Substance HmA[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_selectedForm/N,
              DfH_25degC=DfH_selectedForm/N),
        each amountOfSubstance_start=Kh37*initialHb/(Kh37 + initialH))
        "Deprotonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{92,14},{72,34}})));
      Chemical.Components.Substance HmNH3[
                                         4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=RT*log(Kz25) + DfG_selectedForm/N,
              DfH_25degC=-Hz+DfH_selectedForm/N),
        each amountOfSubstance_start=(initialH^2)*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Protonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{-84,-42},{-64,-22}})));
      Chemical.Components.Reaction z[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{-54,-42},{-34,-22}})));
      Chemical.Components.Substance HmNH2[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_selectedForm/N,
              DfH_25degC=DfH_selectedForm/N),
        each amountOfSubstance_start=initialH*Kz37*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Deprotonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{12,-44},{-8,-24}})));
      Chemical.Components.Reaction c[4](
        each nP=2,
        each KC=KC,
        each nS=2) annotation (Placement(transformation(extent={{20,-42},{40,-22}})));
      Chemical.Interfaces.SubstancePort_b CO2 annotation (Placement(transformation(
              extent={{-8,-26},{12,-6}}), iconTransformation(extent={{10,70},{
                30,90}})));
      Chemical.Components.Substance HmNHCOO[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each amountOfSubstance_start=(Kz37*Kc37*initialCO2)*initialHb/(initialH^2 +
            initialH*Kz37 + Kz37*Kc37*initialCO2),
        each substanceData(DfG_25degC_1bar=DfG_CO2 - RT*log(Kc25) + DfG_selectedForm/
              N,
              DfH_25degC=DfH_CO2+Hc+DfH_selectedForm/N))
        "Carboxylated c site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{70,-44},{50,-24}})));
      Interfaces.SubstancePort_b H annotation (Placement(transformation(extent={{50,-12},
                {70,8}}),      iconTransformation(extent={{70,70},{90,90}})));
    equation

      connect(OxyHm.port_a, o.substrates[1])
        annotation (Line(points={{-68,24},{-60,24},{-60,24}}, color={158,66,200}));
      connect(o.products[1], DeoxyHm.port_a)
        annotation (Line(points={{-40,22},{-28,22},{-28,22}},
                                                     color={158,66,200}));

      for i in 1:4 loop
        connect(h[i].products[2], H) annotation (Line(
            points={{56,20},{60,20},{60,-2}},
            color={158,66,200}));
        connect(speciation.subunitSolution, HmA[i].solution) annotation (Line(
          points={{-4,-56},{-4,-44},{94,-44},{94,12},{88,12},{88,14}},
          color={127,127,0}));
        connect(speciation.subunitSolution, HmAH[i].solution) annotation (Line(
          points={{-4,-56},{-4,-44},{94,-44},{94,12},{14,12}},
          color={127,127,0}));
        connect(HmA[i].port_a, speciation.subunits[i+4]) annotation (Line(
          points={{72,24},{72,-52},{-11,-52},{-11,-51.8}},
          color={158,66,200}));

        connect(o[i].products[2], O2) annotation (Line(points={{-40,26},{-16,26},
                {-16,42},{-18,42}},
                               color={158,66,200}));
        connect(speciation.subunitSolution, DeoxyHm[i].solution) annotation (Line(
              points={{-4,-56},{-4,-44},{94,-44},{94,12},{-12,12}},
                                                              color={127,127,0}));
        connect(speciation.subunitSolution, OxyHm[i].solution) annotation (Line(
              points={{-4,-56},{-4,-44},{94,-44},{94,12},{92,12},{-84,12},{-84,
                14}},                                           color={127,127,
                0}));
        connect(DeoxyHm[i].port_a, speciation.subunits[i]) annotation (Line(
          points={{-28,22},{-28,22},{-12,22},{-12,-22},{-11,-22},{-11,-51.8}},
          color={158,66,200}));

        connect(z[i].products[2], H) annotation (Line(
            points={{-34,-30},{-22,-30},{-22,-2},{60,-2}},
            color={158,66,200}));
        connect(speciation.subunitSolution, HmNH2[i].solution) annotation (Line(
          points={{-4,-56},{-4,-44},{8,-44}},
          color={127,127,0}));
        connect(HmNH2[i].port_a, speciation.subunits[i + 8]) annotation (Line(
          points={{-8,-34},{-11,-34},{-11,-51.8}},
          color={158,66,200}));
        connect(HmNH3[i].solution, speciation.subunitSolution) annotation (Line(
          points={{-80,-42},{-80,-44},{-4,-44},{-4,-56}},
          color={127,127,0}));

        connect(c[i].products[2], H) annotation (Line(
            points={{40,-30},{46,-30},{46,-2},{60,-2}},
            color={158,66,200}));
        connect(CO2, c[i].substrates[2]) annotation (Line(
            points={{2,-16},{16,-16},{16,-30},{20,-30}},
            color={158,66,200}));
        connect(HmNHCOO[i].solution, speciation.subunitSolution) annotation (Line(
          points={{66,-44},{-4,-44},{-4,-56}},
          color={127,127,0}));

      end for;

      connect(speciation.solution, solution) annotation (Line(
          points={{-14,-72},{-22,-72},{-22,-56},{-28,-56},{-28,-72},{-41,-72}},
          color={127,127,0}));
      connect(speciation.port_a, selectedForm) annotation (Line(
          points={{2,-72},{12,-72},{12,-56},{20,-56},{20,-72},{36,-72}},
          color={158,66,200}));
      connect(HmAH.port_a,h. substrates[1]) annotation (Line(
          points={{30,22},{36,22}},
          color={158,66,200}));
      connect(h.products[1],HmA. port_a) annotation (Line(
          points={{56,24},{64,24},{72,24}},
          color={158,66,200}));

      connect(z.products[1], HmNH2.port_a) annotation (Line(
          points={{-34,-34},{-8,-34}},
          color={107,45,134}));

      connect(HmNH3.port_a, z.substrates[1]) annotation (Line(
          points={{-64,-32},{-62,-32},{-60,-32},{-54,-32}},
          color={158,66,200}));

      connect(HmNH2.port_a, c.substrates[1]) annotation (Line(
          points={{-8,-34},{20,-34}},
          color={158,66,200}));

      connect(HmNHCOO.port_a, c.products[1]) annotation (Line(
          points={{50,-34},{40,-34}},
          color={158,66,200}));

      connect(solution, solution) annotation (Line(
          points={{-41,-72},{-41,-72}},
          color={127,127,0}));
      connect(H, H) annotation (Line(
          points={{60,-2},{60,-2}},
          color={158,66,200}));

      annotation (          experiment(StopTime=15000),
        Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<p>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}),
                graphics),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics));
    end HemoglobinQuaternaryForm;

    model HemoglobinMultipleAllostery
      "Multiple-ligand allosteric hemoglobin model"
      extends Modelica.Icons.Example;

      constant Modelica.SIunits.AmountOfSubstance THb = 0.001
        "Total amount of hemoglobin";

      constant Real RT=Modelica.Constants.R*298.15;

      constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        "Amount of solution used for molarity to mole fraction conversion";

      constant Modelica.SIunits.Volume OneLiter = 0.001;

      parameter Real L_old=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.SIunits.Concentration KR=0.000671946
        "oxygen dissociation on relaxed(R) hemoglobin subunit";

      parameter Modelica.SIunits.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.SIunits.MoleFraction KRo37 = KR*OneLiter/AmountOfSolutionIn1L;
      parameter Modelica.SIunits.MoleFraction KTo37 = KT*OneLiter/AmountOfSolutionIn1L;

      parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508);
      parameter Modelica.SIunits.ChemicalPotential DfG_CO2 = -RT*log(0.034/55.508)  + Chemical.Examples.Substances.CarbonDioxide_gas.DfG_25degC_1bar;

      parameter Modelica.SIunits.ChemicalPotential DfG_tT = 0;
      parameter Modelica.SIunits.ChemicalPotential DfG_tR = DfG_tT + RT * log(L);

      parameter Real KC = 1e-3 "Slow down factor";

      parameter Modelica.SIunits.MoleFraction initialO2=1.9594e-07
        "Initial O2 at 37degC, pO2=100Pa";             //at 25degC: 2.342e-8;
      parameter Modelica.SIunits.MoleFraction initialH=10^(-7.2);
      parameter Modelica.SIunits.MoleFraction initialCO2=2.4217e-10
        "Initial CO2 at 37degC, pCO2=40mmHg";          //at 25degC: 3.267e-5;
      parameter Modelica.SIunits.MoleFraction initialCO=1e-10
        "Initial CO at 37degC, pCO=0mmHg";          //at 25degC: 3.267e-5;

      parameter Modelica.SIunits.MoleFraction KRh37 = 10^(-6.89);
      parameter Modelica.SIunits.MoleFraction KTh37 = 10^(-7.52);

      parameter Modelica.SIunits.MoleFraction KRz37 = 10^(-7.25);
      parameter Modelica.SIunits.MoleFraction KTz37 = 10^(-7.73);

      parameter Modelica.SIunits.MoleFraction KRc37 = (10^(-8.35)) / (OneLiter/AmountOfSolutionIn1L);
      parameter Modelica.SIunits.MoleFraction KTc37 = (10^(-7.54)) / (OneLiter/AmountOfSolutionIn1L);

      parameter Real L=L_old
        *
        (((KTh37/((10^(-7.2))+KTh37)) / (KRh37/((10^(-7.2))+KRh37)))^4)
        *
        (((KTz37*((10^(-7.2))^2 + KRz37*(10^(-7.2)) + KRz37*KRc37*(2.4217e-5)))/(KRz37*((10^(-7.2))^2 + KTz37*(10^(-7.2)) + KTz37*KTc37*(2.4217e-5))))^4)
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";

      Chemical.Components.Solution solution(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-56},{100,32}})));

      Chemical.Components.Reaction quaternaryForm(KC=KC)
        annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

      Chemical.Components.Substance O2_free(amountOfSubstance_start=initialO2*
            AmountOfSolutionIn1L, substanceData(DfG_25degC_1bar=DfG_O2,
            DfH_25degC=-11700))
        annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
      Modelica.Blocks.Sources.Clock oxygenSource(offset=1000)
        annotation (Placement(transformation(extent={{-78,48},{-58,68}})));
      Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
          usePartialPressureInput=true, substanceData=Chemical.Examples.Substances.Oxygen_gas,
        Temperature=310.15)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            origin={-36,58})));
      Chemical.Components.GasSolubility partialPressure1(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={-14,32})));

      Chemical.Components.Substance H2O(substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=AmountOfSolutionIn1L - THb - (initialO2 +
            initialCO2)*AmountOfSolutionIn1L)
        annotation (Placement(transformation(extent={{56,-54},{76,-34}})));
      HemoglobinQuaternaryForm                              relaxed(
        Ko37=KRo37,
        DfG_selectedForm=DfG_tR,
        initialO2=initialO2,
        initialHb=THb/(L + 1),
        initialH=initialH,
        Kh37=KRh37,
        Kz37=KRz37,
        Kc37=KRc37,
        initialCO2=initialCO2,
        DfG_O2=DfG_O2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = -41000,
        Hz=8000,
        Hh=127000,
        Kco37=KRo37)
        annotation (Placement(transformation(extent={{-54,-44},{-34,-24}})));
      HemoglobinQuaternaryForm                              tensed(
        Ko37=KTo37,
        DfG_selectedForm=DfG_tT,
        initialO2=initialO2,
        initialHb=THb*L/(L + 1),
        initialH=initialH,
        Kh37=KTh37,
        Kz37=KTz37,
        Kc37=KTc37,
        initialCO2=initialCO2,
        DfG_O2=DfG_O2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = 59000,
        Hz=-51000,
        Hh=59000,
        Kco37=KTo37)
                  annotation (Placement(transformation(extent={{32,-44},{12,-24}})));
      Chemical.Sources.ExternalMoleFraction H(substanceData=Chemical.Examples.Substances.Proton_aqueous,
          MoleFraction=initialH,
        Temperature=310.15)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}},
            origin={-12,-18})));
      Chemical.Components.Substance CO2_free(amountOfSubstance_start=initialCO2
            *AmountOfSolutionIn1L, substanceData(DfG_25degC_1bar=DfG_CO2,
            DfH_25degC=-412900))
        annotation (Placement(transformation(extent={{86,-8},{66,12}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
        PartialPressure(displayUnit="kPa") = 5330,
        Temperature=310.15)                                   annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={62,60})));
      Chemical.Components.GasSolubility partialPressure2(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={62,32})));

      Real sO2 "Hemoglobin oxygen saturation";
      Real sCO2 "Hemoglobin carbon dioxide saturation";
      Real dH "Hemoglobin charge change caused by binding of Bohr's protons";
    equation
      sO2 = (sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x)) /
      (sum(relaxed.DeoxyHm.x) + sum(tensed.DeoxyHm.x) + sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x));

      sCO2 = (sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x)) /
      (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) + sum(relaxed.HmNH2.x) + sum(tensed.HmNH2.x) + sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x));

      dH = (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) - sum(relaxed.HmNHCOO.x) - sum(tensed.HmNHCOO.x) - sum(relaxed.HmA.x) - sum(tensed.HmA.x)) /
      THb;

      connect(oxygenSource.y, oxygen_in_air.partialPressure)
        annotation (Line(points={{-57,58},{-46,58}}, color={0,0,127}));

      connect(oxygen_in_air.port_a, partialPressure1.gas_port) annotation (Line(
          points={{-26,58},{-26,58},{-14,58},{-14,42}},
          color={158,66,200}));
      connect(partialPressure1.liquid_port, O2_free.port_a) annotation (Line(
            points={{-14,22},{-14,-2},{-56,-2}}, color={158,66,200}));

      connect(O2_free.solution, solution.solution) annotation (Line(points={{-72,-12},
              {-72,-54},{60,-54},{60,-55.12}},        color={127,127,0}));
      connect(solution.solution, H2O.solution) annotation (Line(
          points={{60,-55.12},{60,-54}},
          color={127,127,0}));

      connect(relaxed.solution, solution.solution) annotation (Line(
          points={{-48,-42},{-48,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(relaxed.O2, O2_free.port_a) annotation (Line(
          points={{-52,-26},{-52,-2},{-56,-2}},
          color={158,66,200}));
      connect(relaxed.selectedForm, quaternaryForm.substrates[1]) annotation (Line(
          points={{-40,-42},{-22,-42}},
          color={158,66,200}));
      connect(tensed.solution, solution.solution) annotation (Line(
          points={{26,-42},{26,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(tensed.O2, O2_free.port_a) annotation (Line(
          points={{30,-26},{30,-2},{-56,-2}},
          color={158,66,200}));
      connect(tensed.selectedForm, quaternaryForm.products[1]) annotation (Line(
          points={{18,-42},{-2,-42}},
          color={158,66,200}));
      connect(H.port_a, relaxed.H) annotation (Line(
          points={{-22,-18},{-36,-18},{-36,-26}},
          color={158,66,200}));
      connect(H.port_a, tensed.H) annotation (Line(
          points={{-22,-18},{14,-18},{14,-26}},
          color={158,66,200}));
      connect(CO2_gas.port_a, partialPressure2.gas_port) annotation (Line(
          points={{62,50},{62,42}},
          color={158,66,200}));
      connect(partialPressure2.liquid_port, CO2_free.port_a)
        annotation (Line(points={{62,22},{62,2},{66,2}}, color={158,66,200}));
      connect(CO2_free.port_a, tensed.CO2) annotation (Line(
          points={{66,2},{20,2},{20,-26}},
          color={158,66,200}));
      connect(CO2_free.port_a, relaxed.CO2) annotation (Line(
          points={{66,2},{-42,2},{-42,-26}},
          color={158,66,200}));
      connect(CO2_free.solution, solution.solution) annotation (Line(
          points={{82,-8},{82,-8},{82,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      annotation (          experiment(StopTime=15000),
        Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<p>Oxygen dissociation curve of hemoglobin.</p>
<p>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/ODC_pH_CO2_Hb.png\"/></p>
<p>J. W. Severinghaus, &quot;Simple, accurate equations for human blood O2 dissociation computations,&quot; Journal of Applied Physiology, vol. 46, pp. 599-602, 1979.</p>
<p><br>pO2 .. partial pressure of oxygen in gas</p>
<p>pCO2 .. partial pressure of carbon dioxide</p>
<p>sO2 .. oxygen saturation of hemoglobin</p>
<p>pH = log10(aH), where aH is mole fraction based activity of hydrogen ions</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/ODC_T_Hb.png\"/></p>
<p>R. B. Reeves, &quot;The effect of temperature on the oxygen equilibrium curve of human blood,&quot; Respiration physiology, vol. 42, pp. 317-328, 1980.</p>
<p><br>T .. temperature</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end HemoglobinMultipleAllostery;

    model HemoglobinQuaternaryFormCO
      "Hemoglobib quaternary form - part of multiple-ligand allosteric hemoglobin model"

      constant Integer N=12
        "Number of distinguished independent sides in quaternary structure";
      constant Real RT=Modelica.Constants.R*298.15;

      parameter Modelica.SIunits.MolarEnthalpy Ho = 59000
        "Enthalpy of deoxygenation";
      parameter Modelica.SIunits.MoleFraction Ko37 "KRx and KTx at 37degC";
      parameter Modelica.SIunits.MoleFraction Ko25 = Ko37*exp((Ho/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRx and KTx at 25degC";

      parameter Modelica.SIunits.MolarEnthalpy Hco = 59000
        "Enthalpy of carbon monoxide dissociation";
      parameter Modelica.SIunits.MoleFraction Kco37
        "Carboxyhemoglobin dissociation at 37degC";
      parameter Modelica.SIunits.MoleFraction Kco25 = Kco37*exp((Hco/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "Carboxyhemoglobin dissociation at 25degC";

      parameter Modelica.SIunits.MolarEnthalpy Hh
        "Enthalpy of deprotonation of h site";
      parameter Modelica.SIunits.MoleFraction Kh37 "KRhx and KThx at 37 degC";
      parameter Modelica.SIunits.MoleFraction Kh25 = Kh37*exp(((Hh)/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRhx and KThx at 25 degC";

      parameter Modelica.SIunits.MolarEnthalpy Hz
        "Enthalpy of deprotonation of -NH3+ terminus";
      parameter Modelica.SIunits.MoleFraction Kz37 "KRzx and KTzx at 37 degC";
      parameter Modelica.SIunits.MoleFraction Kz25 = Kz37*exp(((Hz)/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRzx and KTzx at 25 degC";

      parameter Modelica.SIunits.MolarEnthalpy Hc "Enthalpy of carboxylation";
      parameter Modelica.SIunits.MoleFraction Kc37 "KRcx and KTcx at 37degC";
      parameter Modelica.SIunits.MoleFraction Kc25 = Kc37*exp((Hc/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRcx and KTcx at 25degC";

      parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508)  + Chemical.Examples.Substances.Oxygen_gas.DfG_25degC_1bar;
      parameter Modelica.SIunits.ChemicalPotential DfH_O2 = Chemical.Examples.Substances.Oxygen_aqueous.DfH_25degC;

      parameter Modelica.SIunits.ChemicalPotential DfG_CO = -RT*log(0.00099/55.508)  + Chemical.Examples.Substances.CarbonMonoxide_gas.DfG_25degC_1bar; //==Chemical.Examples.Substances.CarbonMonoxide_aqueous.DfG_25degC_1bar
      parameter Modelica.SIunits.ChemicalPotential DfH_CO = Chemical.Examples.Substances.CarbonMonoxide_aqueous.DfH_25degC;

      parameter Modelica.SIunits.ChemicalPotential DfG_CO2 = -RT*log(0.034/55.508)  + Chemical.Examples.Substances.CarbonDioxide_gas.DfG_25degC_1bar;
      parameter Modelica.SIunits.ChemicalPotential DfH_CO2 = Chemical.Examples.Substances.CarbonDioxide_aqueous.DfH_25degC;

      parameter Modelica.SIunits.ChemicalPotential DfG_selectedForm
        "DfG_tR and DfG_tT";
      parameter Modelica.SIunits.MolarEnthalpy DfH_selectedForm = 0
        "DfH_tR and DfH_tT";

      parameter Real KC = 1e-3 "Slow down factor";
                               //0.000001
      parameter Modelica.SIunits.MoleFraction initialO2
        "Initial mole fraction of unbound oxygen disoluted around hemoglobin";
      parameter Modelica.SIunits.MoleFraction initialCO
        "Initial mole fraction of unbound carbon monoxide disoluted around hemoglobin";
      parameter Modelica.SIunits.MoleFraction initialH
        "Initial mole fraction of H+";
      parameter Modelica.SIunits.MoleFraction initialCO2
        "Initial mole fraction of unbound carbon dioxide disoluted around hemoglobin";
      parameter Modelica.SIunits.AmountOfSubstance initialHb
        "Initial amount of hemoglobin tetramers in this quaternary form";

      Chemical.Components.Speciation speciation(NumberOfSubunits=N)
        annotation (Placement(transformation(extent={{-18,-72},{2,-52}})));
       // AmountOfSubstance_start=4e-11)
       // AmountOfSubstance_start=totalAmountOfHemoglobin)
      Chemical.Components.Substance OxyHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_O2 + RT*log(Ko25) + DfG_selectedForm/N,
              DfH_25degC=DfH_O2 -Ho + DfH_selectedForm/N),
        each amountOfSubstance_start=(initialO2/Ko37)*initialHb/(1 + initialO2/Ko37 +
            initialCO/Kco37)) "Oxygenated subunit"
        annotation (Placement(transformation(extent={{-88,14},{-68,34}})));

      Chemical.Components.Reaction o[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
      Chemical.Components.Substance DeoxyHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_selectedForm/N,
              DfH_25degC=DfH_selectedForm/N),
        each amountOfSubstance_start=initialHb/(1 + initialO2/Ko37 + initialCO/
            Kco37)) "Deoxygenated subunit"
        annotation (Placement(transformation(extent={{-8,12},{-28,32}})));

      Chemical.Interfaces.SolutionPort solution
        annotation (Placement(transformation(extent={{-50,-82},{-32,-62}}),
            iconTransformation(extent={{-50,-90},{-30,-70}})));
      Chemical.Interfaces.SubstancePort_b O2 annotation (Placement(transformation(
              extent={{-28,32},{-8,52}}),  iconTransformation(extent={{-90,70},{-70,
                90}})));
      Chemical.Interfaces.SubstancePort_a selectedForm
        annotation (Placement(transformation(extent={{26,-82},{46,-62}}),
            iconTransformation(extent={{30,-90},{50,-70}})));
      Chemical.Components.Substance HmAH[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=RT*log(Kh25) + DfG_selectedForm/N,
              DfH_25degC=-Hh+DfH_selectedForm/N),
        each amountOfSubstance_start=initialH*initialHb/(Kh37 + initialH))
        "Protonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{10,12},{30,32}})));
      Chemical.Components.Reaction h[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{36,32},{56,12}})));
      Chemical.Components.Substance HmA[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_selectedForm/N,
              DfH_25degC=DfH_selectedForm/N),
        each amountOfSubstance_start=Kh37*initialHb/(Kh37 + initialH))
        "Deprotonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{92,14},{72,34}})));
      Chemical.Components.Substance HmNH3[
                                         4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=RT*log(Kz25) + DfG_selectedForm/N,
              DfH_25degC=-Hz+DfH_selectedForm/N),
        each amountOfSubstance_start=(initialH^2)*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Protonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{-84,-42},{-64,-22}})));
      Chemical.Components.Reaction z[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{-54,-42},{-34,-22}})));
      Chemical.Components.Substance HmNH2[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_selectedForm/N,
              DfH_25degC=DfH_selectedForm/N),
        each amountOfSubstance_start=initialH*Kz37*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Deprotonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{12,-44},{-8,-24}})));
      Chemical.Components.Reaction c[4](
        each nP=2,
        each KC=KC,
        each nS=2) annotation (Placement(transformation(extent={{20,-42},{40,-22}})));
      Chemical.Interfaces.SubstancePort_b CO2 annotation (Placement(transformation(
              extent={{-8,-26},{12,-6}}), iconTransformation(extent={{10,70},{
                30,90}})));
      Chemical.Components.Substance HmNHCOO[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each amountOfSubstance_start=(Kz37*Kc37*initialCO2)*initialHb/(initialH^2 +
            initialH*Kz37 + Kz37*Kc37*initialCO2),
        each substanceData(DfG_25degC_1bar=DfG_CO2 - RT*log(Kc25) + DfG_selectedForm/
              N,
              DfH_25degC=DfH_CO2+Hc+DfH_selectedForm/N))
        "Carboxylated c site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{70,-44},{50,-24}})));
      Interfaces.SubstancePort_b H annotation (Placement(transformation(extent={{50,-12},
                {70,8}}),      iconTransformation(extent={{70,70},{90,90}})));
      Chemical.Components.Substance COHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each amountOfSubstance_start=(initialCO/Kco37)*initialHb/(1 + initialO2/
            Ko37 + initialCO/Kco37),
        each substanceData(DfH_25degC=DfH_CO - Hco + DfH_selectedForm/N,
            DfG_25degC_1bar=DfG_CO + RT*log(Kco25) + DfG_selectedForm/N))
        "Subunit with Carbon Monoxide"
        annotation (Placement(transformation(extent={{78,44},{58,64}})));
      Chemical.Components.Reaction o1[
                                     4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{46,44},{26,64}})));
      Chemical.Interfaces.SubstancePort_b CO annotation (Placement(transformation(
              extent={{-12,62},{8,82}}),   iconTransformation(extent={{-50,70},
                {-30,90}})));
    equation

      connect(OxyHm.port_a, o.substrates[1])
        annotation (Line(points={{-68,24},{-60,24}},          color={158,66,200}));
      connect(o.products[1], DeoxyHm.port_a)
        annotation (Line(points={{-40,22},{-28,22}}, color={158,66,200}));

      for i in 1:4 loop
        connect(h[i].products[2], H) annotation (Line(
            points={{56,20},{60,20},{60,-2}},
            color={158,66,200}));
        connect(speciation.subunitSolution, HmA[i].solution) annotation (Line(
          points={{-4,-56},{-4,-44},{94,-44},{94,12},{88,12},{88,14}},
          color={127,127,0}));
        connect(speciation.subunitSolution, HmAH[i].solution) annotation (Line(
          points={{-4,-56},{-4,-44},{94,-44},{94,12},{14,12}},
          color={127,127,0}));
        connect(HmA[i].port_a, speciation.subunits[i+4]) annotation (Line(
          points={{72,24},{72,-52},{-11,-52},{-11,-51.8}},
          color={158,66,200}));

        connect(o[i].products[2], O2) annotation (Line(points={{-40,26},{-16,26},{-16,
                42},{-18,42}}, color={158,66,200}));
        connect(speciation.subunitSolution, DeoxyHm[i].solution) annotation (Line(
              points={{-4,-56},{-4,-44},{94,-44},{94,12},{-12,12}},
                                                              color={127,127,0}));
        connect(speciation.subunitSolution, OxyHm[i].solution) annotation (Line(
              points={{-4,-56},{-4,-44},{94,-44},{94,12},{92,12},{-84,12},{-84,14}},
                                                                color={127,127,
                0}));
        connect(DeoxyHm[i].port_a, speciation.subunits[i]) annotation (Line(
          points={{-28,22},{-28,22},{-12,22},{-12,-22},{-11,-22},{-11,-51.8}},
          color={158,66,200}));

        connect(z[i].products[2], H) annotation (Line(
            points={{-34,-30},{-22,-30},{-22,-2},{60,-2}},
            color={158,66,200}));
        connect(speciation.subunitSolution, HmNH2[i].solution) annotation (Line(
          points={{-4,-56},{-4,-44},{8,-44}},
          color={127,127,0}));
        connect(HmNH2[i].port_a, speciation.subunits[i + 8]) annotation (Line(
          points={{-8,-34},{-11,-34},{-11,-51.8}},
          color={158,66,200}));
        connect(HmNH3[i].solution, speciation.subunitSolution) annotation (Line(
          points={{-80,-42},{-80,-44},{-4,-44},{-4,-56}},
          color={127,127,0}));

        connect(c[i].products[2], H) annotation (Line(
            points={{40,-30},{46,-30},{46,-2},{60,-2}},
            color={158,66,200}));
        connect(CO2, c[i].substrates[2]) annotation (Line(
            points={{2,-16},{16,-16},{16,-30},{20,-30}},
            color={158,66,200}));
        connect(HmNHCOO[i].solution, speciation.subunitSolution) annotation (Line(
          points={{66,-44},{-4,-44},{-4,-56}},
          color={127,127,0}));

        connect(COHm[i].solution, speciation.subunitSolution) annotation (Line(
          points={{74,44},{74,42},{94,42},{94,-44},{-4,-44},{-4,-56}},
          color={127,127,0},
          smooth=Smooth.None));
        connect(o1[i].products[2], CO) annotation (Line(
          points={{26,56},{14,56},{14,72},{-2,72}},
          color={158,66,200},
          smooth=Smooth.None));

      end for;

      connect(speciation.solution, solution) annotation (Line(
          points={{-14,-72},{-22,-72},{-22,-56},{-28,-56},{-28,-72},{-41,-72}},
          color={127,127,0}));
      connect(speciation.port_a, selectedForm) annotation (Line(
          points={{2,-72},{12,-72},{12,-56},{20,-56},{20,-72},{36,-72}},
          color={158,66,200}));
      connect(HmAH.port_a,h. substrates[1]) annotation (Line(
          points={{30,22},{36,22}},
          color={158,66,200}));
      connect(h.products[1],HmA. port_a) annotation (Line(
          points={{56,24},{64,24},{72,24}},
          color={158,66,200}));

      connect(z.products[1], HmNH2.port_a) annotation (Line(
          points={{-34,-34},{-8,-34}},
          color={107,45,134}));

      connect(HmNH3.port_a, z.substrates[1]) annotation (Line(
          points={{-64,-32},{-62,-32},{-60,-32},{-54,-32}},
          color={158,66,200}));

      connect(HmNH2.port_a, c.substrates[1]) annotation (Line(
          points={{-8,-34},{20,-34}},
          color={158,66,200}));

      connect(HmNHCOO.port_a, c.products[1]) annotation (Line(
          points={{50,-34},{40,-34}},
          color={158,66,200}));

      connect(solution, solution) annotation (Line(
          points={{-41,-72},{-41,-72}},
          color={127,127,0}));
      connect(H, H) annotation (Line(
          points={{60,-2},{60,-2}},
          color={158,66,200}));

      connect(COHm.port_a, o1.substrates[1]) annotation (Line(
          points={{58,54},{46,54}},
          color={158,66,200},
          smooth=Smooth.None));

      connect(DeoxyHm.port_a, o1.products[1]) annotation (Line(
          points={{-28,22},{-2,22},{-2,52},{26,52}},
          color={158,66,200},
          smooth=Smooth.None));
      annotation (          experiment(StopTime=15000),
        Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<p>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}),
                graphics),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics));
    end HemoglobinQuaternaryFormCO;

    model HemoglobinMultipleAllosteryCO
      "Multiple-ligand allosteric hemoglobin model"
      extends Modelica.Icons.Example;

      constant Modelica.SIunits.AmountOfSubstance THb = 0.001
        "Total amount of hemoglobin";

      constant Real RT=Modelica.Constants.R*298.15;

      constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        "Amount of solution used for molarity to mole fraction conversion";

      constant Modelica.SIunits.Volume OneLiter = 0.001;

      parameter Real L_old=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.SIunits.Concentration KR=0.000671946
        "oxygen dissociation on relaxed(R) hemoglobin subunit";

      parameter Modelica.SIunits.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.SIunits.MoleFraction KRo37 = KR*OneLiter/AmountOfSolutionIn1L;
      parameter Modelica.SIunits.MoleFraction KTo37 = KT*OneLiter/AmountOfSolutionIn1L;

      parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508);
      parameter Modelica.SIunits.ChemicalPotential DfG_CO2 = -RT*log(0.034/55.508)  + Chemical.Examples.Substances.CarbonDioxide_gas.DfG_25degC_1bar;

      parameter Modelica.SIunits.ChemicalPotential DfG_tT = 0;
      parameter Modelica.SIunits.ChemicalPotential DfG_tR = DfG_tT + RT * log(L);

      parameter Real KC = 1e-3 "Slow down factor";

      parameter Modelica.SIunits.MoleFraction initialO2=1.9594e-07
        "Initial O2 at 37degC, pO2=100Pa";             //at 25degC: 2.342e-8;
      parameter Modelica.SIunits.MoleFraction initialH=10^(-7.2);
      parameter Modelica.SIunits.MoleFraction initialCO2=2.4217e-10
        "Initial CO2 at 37degC, pCO2=40mmHg";          //at 25degC: 3.267e-5;
      parameter Modelica.SIunits.MoleFraction initialCO=1e-10
        "Initial CO at 37degC, pCO=0mmHg";          //at 25degC: 3.267e-5;

      parameter Modelica.SIunits.MoleFraction KRh37 = 10^(-6.89);
      parameter Modelica.SIunits.MoleFraction KTh37 = 10^(-7.52);

      parameter Modelica.SIunits.MoleFraction KRz37 = 10^(-7.25);
      parameter Modelica.SIunits.MoleFraction KTz37 = 10^(-7.73);

      parameter Modelica.SIunits.MoleFraction KRc37 = (10^(-8.35)) / (OneLiter/AmountOfSolutionIn1L);
      parameter Modelica.SIunits.MoleFraction KTc37 = (10^(-7.54)) / (OneLiter/AmountOfSolutionIn1L);

      parameter Real L=L_old
        *
        (((KTh37/((10^(-7.2))+KTh37)) / (KRh37/((10^(-7.2))+KRh37)))^4)
        *
        (((KTz37*((10^(-7.2))^2 + KRz37*(10^(-7.2)) + KRz37*KRc37*(2.4217e-5)))/(KRz37*((10^(-7.2))^2 + KTz37*(10^(-7.2)) + KTz37*KTc37*(2.4217e-5))))^4)
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";

      Chemical.Components.Solution solution(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-56},{100,32}})));

      Chemical.Components.Reaction quaternaryForm(KC=KC)
        annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

      Chemical.Components.Substance O2_free(amountOfSubstance_start=initialO2*
            AmountOfSolutionIn1L, substanceData(DfG_25degC_1bar=DfG_O2,
            DfH_25degC=-11700))
        annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
      Modelica.Blocks.Sources.Clock oxygenSource(offset=1000)
        annotation (Placement(transformation(extent={{-78,48},{-58,68}})));
      Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
          usePartialPressureInput=true, substanceData=Chemical.Examples.Substances.Oxygen_gas,
        Temperature=310.15)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            origin={-36,58})));
      Chemical.Components.GasSolubility partialPressure1(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={-14,32})));

      Chemical.Components.Substance H2O(substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=AmountOfSolutionIn1L - THb - (initialO2 +
            initialCO2)*AmountOfSolutionIn1L)
        annotation (Placement(transformation(extent={{56,-54},{76,-34}})));
      HemoglobinQuaternaryFormCO                            relaxed(
        DfG_selectedForm=DfG_tR,
        initialO2=initialO2,
        initialHb=THb/(L + 1),
        initialH=initialH,
        Kh37=KRh37,
        Kz37=KRz37,
        Kc37=KRc37,
        initialCO2=initialCO2,
        DfG_O2=DfG_O2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = -41000,
        Hz=8000,
        Hh=127000,
        initialCO=initialCO,
        Ko37=KRo37,
        Kco37=KRo37/3200)
        annotation (Placement(transformation(extent={{-54,-44},{-34,-24}})));
      HemoglobinQuaternaryFormCO                            tensed(
        DfG_selectedForm=DfG_tT,
        initialO2=initialO2,
        initialHb=THb*L/(L + 1),
        initialH=initialH,
        Kh37=KTh37,
        Kz37=KTz37,
        Kc37=KTc37,
        initialCO2=initialCO2,
        DfG_O2=DfG_O2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = 59000,
        Hz=-51000,
        Hh=59000,
        initialCO=initialCO,
        Ko37=KTo37,
        Kco37=KTo37/3200)
                  annotation (Placement(transformation(extent={{32,-44},{12,-24}})));
      Chemical.Sources.ExternalMoleFraction H(substanceData=Chemical.Examples.Substances.Proton_aqueous,
          MoleFraction=initialH,
        Temperature=310.15)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}},
            origin={-12,-18})));
      Chemical.Components.Substance CO2_free(amountOfSubstance_start=initialCO2
            *AmountOfSolutionIn1L, substanceData(DfG_25degC_1bar=DfG_CO2,
            DfH_25degC=-412900))
        annotation (Placement(transformation(extent={{86,-8},{66,12}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
        PartialPressure(displayUnit="kPa") = 5330,
        Temperature=310.15)                                   annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={62,60})));
      Chemical.Components.GasSolubility partialPressure2(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={62,32})));

      Real sCO "Hemoglobin carbon monoxide saturation";
      Real sO2 "Hemoglobin oxygen saturation";
      Real sCO2 "Hemoglobin carbon dioxide saturation";
      Real dH "Hemoglobin charge change caused by binding of Bohr's protons";
      Chemical.Components.Substance CO_free(substanceData=Chemical.Examples.Substances.CarbonMonoxide_aqueous,
          amountOfSubstance_start=initialCO*AmountOfSolutionIn1L)
        annotation (Placement(transformation(extent={{-92,8},{-72,28}})));
      Chemical.Components.GasSolubility partialPressure3(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={26,32})));
      Chemical.Sources.ExternalIdealGasSubstance CO_gas(
        substanceData=Chemical.Examples.Substances.CarbonMonoxide_gas,
        PartialPressure(displayUnit="Pa") = 1e-3,
        Temperature=310.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={26,60})));
    equation
      sCO = (sum(relaxed.COHm.x) + sum(tensed.COHm.x)) /
      (sum(relaxed.DeoxyHm.x) + sum(tensed.DeoxyHm.x) + sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x)+ sum(relaxed.COHm.x) + sum(tensed.COHm.x));

      sO2 = (sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x)) /
      (sum(relaxed.DeoxyHm.x) + sum(tensed.DeoxyHm.x) + sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x)+ sum(relaxed.COHm.x) + sum(tensed.COHm.x));

      sCO2 = (sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x)) /
      (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) + sum(relaxed.HmNH2.x) + sum(tensed.HmNH2.x) + sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x));

      dH = (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) - sum(relaxed.HmNHCOO.x) - sum(tensed.HmNHCOO.x) - sum(relaxed.HmA.x) - sum(tensed.HmA.x)) /
      THb;

      connect(oxygenSource.y, oxygen_in_air.partialPressure)
        annotation (Line(points={{-57,58},{-46,58}}, color={0,0,127}));

      connect(oxygen_in_air.port_a, partialPressure1.gas_port) annotation (Line(
          points={{-26,58},{-26,58},{-14,58},{-14,42}},
          color={158,66,200}));
      connect(partialPressure1.liquid_port, O2_free.port_a) annotation (Line(
            points={{-14,22},{-14,-2},{-56,-2}}, color={158,66,200}));

      connect(O2_free.solution, solution.solution) annotation (Line(points={{-72,-12},
              {-72,-54},{60,-54},{60,-55.12}},        color={127,127,0}));
      connect(solution.solution, H2O.solution) annotation (Line(
          points={{60,-55.12},{60,-54}},
          color={127,127,0}));

      connect(relaxed.solution, solution.solution) annotation (Line(
          points={{-48,-42},{-48,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(relaxed.O2, O2_free.port_a) annotation (Line(
          points={{-52,-26},{-52,-2},{-56,-2}},
          color={158,66,200}));
      connect(relaxed.selectedForm, quaternaryForm.substrates[1]) annotation (Line(
          points={{-40,-42},{-22,-42}},
          color={158,66,200}));
      connect(tensed.solution, solution.solution) annotation (Line(
          points={{26,-42},{26,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(tensed.O2, O2_free.port_a) annotation (Line(
          points={{30,-26},{30,-2},{-56,-2}},
          color={158,66,200}));
      connect(tensed.selectedForm, quaternaryForm.products[1]) annotation (Line(
          points={{18,-42},{-2,-42}},
          color={158,66,200}));
      connect(H.port_a, relaxed.H) annotation (Line(
          points={{-22,-18},{-36,-18},{-36,-26}},
          color={158,66,200}));
      connect(H.port_a, tensed.H) annotation (Line(
          points={{-22,-18},{14,-18},{14,-26}},
          color={158,66,200}));
      connect(CO2_gas.port_a, partialPressure2.gas_port) annotation (Line(
          points={{62,50},{62,42}},
          color={158,66,200}));
      connect(partialPressure2.liquid_port, CO2_free.port_a)
        annotation (Line(points={{62,22},{62,2},{66,2}}, color={158,66,200}));
      connect(CO2_free.port_a, tensed.CO2) annotation (Line(
          points={{66,2},{20,2},{20,-26}},
          color={158,66,200}));
      connect(CO2_free.port_a, relaxed.CO2) annotation (Line(
          points={{66,2},{-42,2},{-42,-26}},
          color={158,66,200}));
      connect(CO2_free.solution, solution.solution) annotation (Line(
          points={{82,-8},{82,-8},{82,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(CO_free.solution, solution.solution) annotation (Line(
          points={{-88,8},{-88,-14},{-72,-14},{-72,-54},{60,-54},{60,-55.12}},
          color={127,127,0},
          smooth=Smooth.None));
      connect(CO_gas.port_a, partialPressure3.gas_port) annotation (Line(
          points={{26,50},{26,42}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(partialPressure3.liquid_port, CO_free.port_a) annotation (Line(
          points={{26,22},{26,18},{-72,18}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(CO_free.port_a, relaxed.CO) annotation (Line(
          points={{-72,18},{-48,18},{-48,-26}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(CO_free.port_a, tensed.CO) annotation (Line(
          points={{-72,18},{26,18},{26,-26}},
          color={158,66,200},
          smooth=Smooth.None));
      annotation (          experiment(StopTime=15000),
        Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<p>Oxygen dissociation curve of hemoglobin.</p>
<p>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/ODC_pH_CO2_Hb.png\"/></p>
<p>J. W. Severinghaus, &quot;Simple, accurate equations for human blood O2 dissociation computations,&quot; Journal of Applied Physiology, vol. 46, pp. 599-602, 1979.</p>
<p><br>pO2 .. partial pressure of oxygen in gas</p>
<p>pCO2 .. partial pressure of carbon dioxide</p>
<p>sO2 .. oxygen saturation of hemoglobin</p>
<p>pH = log10(aH), where aH is mole fraction based activity of hydrogen ions</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/ODC_T_Hb.png\"/></p>
<p>R. B. Reeves, &quot;The effect of temperature on the oxygen equilibrium curve of human blood,&quot; Respiration physiology, vol. 42, pp. 317-328, 1980.</p>
<p><br>T .. temperature</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end HemoglobinMultipleAllosteryCO;

    model HemoglobinQuaternaryForm_
      "Hemoglobib quaternary form - part of multiple-ligand allosteric hemoglobin model"

      constant Integer N=12
        "Number of distinguished independent sides in quaternary structure";
      constant Real RT=Modelica.Constants.R*298.15;
    /*
  parameter Modelica.SIunits.MolarEnthalpy Ho = 59000 
    "Enthalpy of deoxygenation";
  parameter Modelica.SIunits.MoleFraction Ko37 "KRx and KTx at 37degC";
  parameter Modelica.SIunits.MoleFraction Ko25 = Ko37*exp((Ho/Modelica.Constants.R)*(1/310.15 - 1/298.15)) 
    "KRx and KTx at 25degC";
*/
      parameter Modelica.SIunits.MolarEnthalpy Hco = 59000
        "Enthalpy of carbon monoxide dissociation";
      parameter Modelica.SIunits.MoleFraction Kco37
        "Carboxyhemoglobin dissociation at 37degC";
      parameter Modelica.SIunits.MoleFraction Kco25 = Kco37*exp((Hco/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "Carboxyhemoglobin dissociation at 25degC";

      parameter Modelica.SIunits.MolarEnthalpy Hh
        "Enthalpy of deprotonation of h site";
      parameter Modelica.SIunits.MoleFraction Kh37 "KRhx and KThx at 37 degC";
      parameter Modelica.SIunits.MoleFraction Kh25 = Kh37*exp(((Hh)/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRhx and KThx at 25 degC";

      parameter Modelica.SIunits.MolarEnthalpy Hz
        "Enthalpy of deprotonation of -NH3+ terminus";
      parameter Modelica.SIunits.MoleFraction Kz37 "KRzx and KTzx at 37 degC";
      parameter Modelica.SIunits.MoleFraction Kz25 = Kz37*exp(((Hz)/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRzx and KTzx at 25 degC";

      parameter Modelica.SIunits.MolarEnthalpy Hc "Enthalpy of carboxylation";
      parameter Modelica.SIunits.MoleFraction Kc37 "KRcx and KTcx at 37degC";
      parameter Modelica.SIunits.MoleFraction Kc25 = Kc37*exp((Hc/Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRcx and KTcx at 25degC";

    //  parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508)  + Chemical.Examples.Substances.Oxygen_gas.DfG_25degC_1bar;
    //  parameter Modelica.SIunits.ChemicalPotential DfH_O2 = Chemical.Examples.Substances.Oxygen_aqueous.DfH_25degC;

      parameter Modelica.SIunits.ChemicalPotential DfG_CO = -RT*log(0.00099/55.508)  + Chemical.Examples.Substances.CarbonMonoxide_gas.DfG_25degC_1bar; //==Chemical.Examples.Substances.CarbonMonoxide_aqueous.DfG_25degC_1bar
      parameter Modelica.SIunits.ChemicalPotential DfH_CO = Chemical.Examples.Substances.CarbonMonoxide_aqueous.DfH_25degC;

      parameter Modelica.SIunits.ChemicalPotential DfG_CO2 = -RT*log(0.034/55.508)  + Chemical.Examples.Substances.CarbonDioxide_gas.DfG_25degC_1bar;
      parameter Modelica.SIunits.ChemicalPotential DfH_CO2 = Chemical.Examples.Substances.CarbonDioxide_aqueous.DfH_25degC;

      parameter Modelica.SIunits.ChemicalPotential DfG_selectedForm
        "DfG_tR and DfG_tT";
      parameter Modelica.SIunits.MolarEnthalpy DfH_selectedForm = 0
        "DfH_tR and DfH_tT";

      parameter Real KC = 1e-3 "Slow down factor";
                               //0.000001
    //  parameter Modelica.SIunits.MoleFraction initialO2
    //    "Initial mole fraction of unbound oxygen disoluted around hemoglobin";
      parameter Modelica.SIunits.MoleFraction initialCO
        "Initial mole fraction of unbound oxygen disoluted around hemoglobin";
      parameter Modelica.SIunits.MoleFraction initialH
        "Initial mole fraction of H+";
      parameter Modelica.SIunits.MoleFraction initialCO2
        "Initial mole fraction of unbound carbon dioxide disoluted around hemoglobin";
      parameter Modelica.SIunits.AmountOfSubstance initialHb
        "Initial amount of hemoglobin tetramers in this quaternary form";

      Chemical.Components.Speciation speciation(NumberOfSubunits=N)
        annotation (Placement(transformation(extent={{-18,-72},{2,-52}})));
       // AmountOfSubstance_start=4e-11)
       // AmountOfSubstance_start=totalAmountOfHemoglobin)
      Chemical.Components.Substance COHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfH_25degC=DfH_CO - Hco + DfH_selectedForm/N,
            DfG_25degC_1bar=DfG_CO + RT*log(Kco25) + DfG_selectedForm/N),
        each amountOfSubstance_start=initialCO*initialHb/(Kco37 + initialCO))
        "Carboxysubunit"
        annotation (Placement(transformation(extent={{-88,14},{-68,34}})));

      Chemical.Components.Reaction o[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
      Chemical.Components.Substance DeoxyHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_selectedForm/N,
              DfH_25degC=DfH_selectedForm/N),
        each amountOfSubstance_start=Kco37*initialHb/(Kco37 + initialCO))
        "Deoxygenated subunit"
        annotation (Placement(transformation(extent={{-8,12},{-28,32}})));

      Chemical.Interfaces.SolutionPort solution
        annotation (Placement(transformation(extent={{-50,-82},{-32,-62}}),
            iconTransformation(extent={{-50,-90},{-30,-70}})));
      Chemical.Interfaces.SubstancePort_b CO annotation (Placement(transformation(
              extent={{-28,32},{-8,52}}),  iconTransformation(extent={{-90,70},{-70,
                90}})));
      Chemical.Interfaces.SubstancePort_a selectedForm
        annotation (Placement(transformation(extent={{26,-82},{46,-62}}),
            iconTransformation(extent={{30,-90},{50,-70}})));
      Chemical.Components.Substance HmAH[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=RT*log(Kh25) + DfG_selectedForm/N,
              DfH_25degC=-Hh+DfH_selectedForm/N),
        each amountOfSubstance_start=initialH*initialHb/(Kh37 + initialH))
        "Protonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{10,12},{30,32}})));
      Chemical.Components.Reaction h[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{36,32},{56,12}})));
      Chemical.Components.Substance HmA[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_selectedForm/N,
              DfH_25degC=DfH_selectedForm/N),
        each amountOfSubstance_start=Kh37*initialHb/(Kh37 + initialH))
        "Deprotonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{92,14},{72,34}})));
      Chemical.Components.Substance HmNH3[
                                         4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=RT*log(Kz25) + DfG_selectedForm/N,
              DfH_25degC=-Hz+DfH_selectedForm/N),
        each amountOfSubstance_start=(initialH^2)*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Protonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{-84,-42},{-64,-22}})));
      Chemical.Components.Reaction z[4](each nP=2, each KC=KC)
        annotation (Placement(transformation(extent={{-54,-42},{-34,-22}})));
      Chemical.Components.Substance HmNH2[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG_25degC_1bar=DfG_selectedForm/N,
              DfH_25degC=DfH_selectedForm/N),
        each amountOfSubstance_start=initialH*Kz37*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Deprotonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{12,-44},{-8,-24}})));
      Chemical.Components.Reaction c[4](
        each nP=2,
        each KC=KC,
        each nS=2) annotation (Placement(transformation(extent={{20,-42},{40,-22}})));
      Chemical.Interfaces.SubstancePort_b CO2 annotation (Placement(transformation(
              extent={{-8,-26},{12,-6}}), iconTransformation(extent={{10,70},{
                30,90}})));
      Chemical.Components.Substance HmNHCOO[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each amountOfSubstance_start=(Kz37*Kc37*initialCO2)*initialHb/(initialH^2 +
            initialH*Kz37 + Kz37*Kc37*initialCO2),
        each substanceData(DfG_25degC_1bar=DfG_CO2 - RT*log(Kc25) + DfG_selectedForm/
              N,
              DfH_25degC=DfH_CO2+Hc+DfH_selectedForm/N))
        "Carboxylated c site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{70,-44},{50,-24}})));
      Interfaces.SubstancePort_b H annotation (Placement(transformation(extent={{50,-12},
                {70,8}}),      iconTransformation(extent={{70,70},{90,90}})));
    equation

      connect(COHm.port_a, o.substrates[1])
        annotation (Line(points={{-68,24},{-60,24}}, color={158,66,200}));
      connect(o.products[1], DeoxyHm.port_a)
        annotation (Line(points={{-40,22},{-28,22},{-28,22}},
                                                     color={158,66,200}));

      for i in 1:4 loop
        connect(h[i].products[2], H) annotation (Line(
            points={{56,20},{60,20},{60,-2}},
            color={158,66,200}));
        connect(speciation.subunitSolution, HmA[i].solution) annotation (Line(
          points={{-4,-56},{-4,-44},{94,-44},{94,12},{88,12},{88,14}},
          color={127,127,0}));
        connect(speciation.subunitSolution, HmAH[i].solution) annotation (Line(
          points={{-4,-56},{-4,-44},{94,-44},{94,12},{14,12}},
          color={127,127,0}));
        connect(HmA[i].port_a, speciation.subunits[i+4]) annotation (Line(
          points={{72,24},{72,-52},{-11,-52},{-11,-51.8}},
          color={158,66,200}));

        connect(o[i].products[2],CO)  annotation (Line(points={{-40,26},{-32,26},{-32,
                42},{-18,42}}, color={158,66,200}));
        connect(speciation.subunitSolution, DeoxyHm[i].solution) annotation (Line(
              points={{-4,-56},{-4,-44},{94,-44},{94,12},{-12,12}},
                                                              color={127,127,0}));
        connect(speciation.subunitSolution, COHm[i].solution) annotation (Line(
              points={{-4,-56},{-4,-44},{94,-44},{94,12},{92,12},{-84,12},{-84,14}},
              color={127,127,0}));
        connect(DeoxyHm[i].port_a, speciation.subunits[i]) annotation (Line(
          points={{-28,22},{-28,22},{-12,22},{-12,-22},{-11,-22},{-11,-51.8}},
          color={158,66,200}));

        connect(z[i].products[2], H) annotation (Line(
            points={{-34,-30},{-22,-30},{-22,-2},{60,-2}},
            color={158,66,200}));
        connect(speciation.subunitSolution, HmNH2[i].solution) annotation (Line(
          points={{-4,-56},{-4,-44},{8,-44}},
          color={127,127,0}));
        connect(HmNH2[i].port_a, speciation.subunits[i + 8]) annotation (Line(
          points={{-8,-34},{-11,-34},{-11,-51.8}},
          color={158,66,200}));
        connect(HmNH3[i].solution, speciation.subunitSolution) annotation (Line(
          points={{-80,-42},{-80,-44},{-4,-44},{-4,-56}},
          color={127,127,0}));

        connect(c[i].products[2], H) annotation (Line(
            points={{40,-30},{46,-30},{46,-2},{60,-2}},
            color={158,66,200}));
        connect(CO2, c[i].substrates[2]) annotation (Line(
            points={{2,-16},{16,-16},{16,-30},{20,-30}},
            color={158,66,200}));
        connect(HmNHCOO[i].solution, speciation.subunitSolution) annotation (Line(
          points={{66,-44},{-4,-44},{-4,-56}},
          color={127,127,0}));

      end for;

      connect(speciation.solution, solution) annotation (Line(
          points={{-14,-72},{-22,-72},{-22,-56},{-28,-56},{-28,-72},{-41,-72}},
          color={127,127,0}));
      connect(speciation.port_a, selectedForm) annotation (Line(
          points={{2,-72},{12,-72},{12,-56},{20,-56},{20,-72},{36,-72}},
          color={158,66,200}));
      connect(HmAH.port_a,h. substrates[1]) annotation (Line(
          points={{30,22},{36,22}},
          color={158,66,200}));
      connect(h.products[1],HmA. port_a) annotation (Line(
          points={{56,24},{64,24},{72,24}},
          color={158,66,200}));

      connect(z.products[1], HmNH2.port_a) annotation (Line(
          points={{-34,-34},{-8,-34}},
          color={107,45,134}));

      connect(HmNH3.port_a, z.substrates[1]) annotation (Line(
          points={{-64,-32},{-62,-32},{-60,-32},{-54,-32}},
          color={158,66,200}));

      connect(HmNH2.port_a, c.substrates[1]) annotation (Line(
          points={{-8,-34},{20,-34}},
          color={158,66,200}));

      connect(HmNHCOO.port_a, c.products[1]) annotation (Line(
          points={{50,-34},{40,-34}},
          color={158,66,200}));

      connect(solution, solution) annotation (Line(
          points={{-41,-72},{-41,-72}},
          color={127,127,0}));
      connect(H, H) annotation (Line(
          points={{60,-2},{60,-2}},
          color={158,66,200}));

      annotation (          experiment(StopTime=15000),
        Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<p>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}),
                graphics),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics));
    end HemoglobinQuaternaryForm_;

    model HemoglobinMultipleAllostery_
      "Multiple-ligand allosteric hemoglobin model"
      extends Modelica.Icons.Example;

      constant Modelica.SIunits.AmountOfSubstance THb = 0.001
        "Total amount of hemoglobin";

      constant Real RT=Modelica.Constants.R*298.15;

      constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        "Amount of solution used for molarity to mole fraction conversion";

      constant Modelica.SIunits.Volume OneLiter = 0.001;

      parameter Real KRco37=KRo37/3200;
      parameter Real KTco37=KTo37/3200;

      parameter Real pH=7.4;
      parameter Modelica.SIunits.Pressure pCO2=5330;

      parameter Real L_old=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.SIunits.Concentration KR=0.000671946
        "oxygen dissociation on relaxed(R) hemoglobin subunit";

      parameter Modelica.SIunits.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.SIunits.MoleFraction KRo37 = KR*OneLiter/AmountOfSolutionIn1L;
      parameter Modelica.SIunits.MoleFraction KTo37 = KT*OneLiter/AmountOfSolutionIn1L;

    //  parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508);
      parameter Modelica.SIunits.ChemicalPotential DfG_CO2 = -RT*log(0.034/55.508)  + Chemical.Examples.Substances.CarbonDioxide_gas.DfG_25degC_1bar;

      parameter Modelica.SIunits.ChemicalPotential DfG_tT = 0;
      parameter Modelica.SIunits.ChemicalPotential DfG_tR = DfG_tT + RT * log(L);

      parameter Real KC = 1e-3 "Slow down factor";

      parameter Modelica.SIunits.MoleFraction initialO2=1.9594e-07
        "Initial O2 at 37degC, pO2=100Pa";             //at 25degC: 2.342e-8;
      parameter Modelica.SIunits.MoleFraction initialH=10^(-pH);
      parameter Modelica.SIunits.MoleFraction initialCO2=2.4217e-10
        "Initial CO2 at 37degC, pCO2=40mmHg";          //at 25degC: 3.267e-5;
      parameter Modelica.SIunits.MoleFraction initialCO=1e-10
        "Initial CO at 37degC, pCO=0mmHg";          //at 25degC: 3.267e-5;

      parameter Modelica.SIunits.MoleFraction KRh37 = 10^(-6.89);
      parameter Modelica.SIunits.MoleFraction KTh37 = 10^(-7.52);

      parameter Modelica.SIunits.MoleFraction KRz37 = 10^(-7.25);
      parameter Modelica.SIunits.MoleFraction KTz37 = 10^(-7.73);

      parameter Modelica.SIunits.MoleFraction KRc37 = (10^(-8.35)) / (OneLiter/AmountOfSolutionIn1L);
      parameter Modelica.SIunits.MoleFraction KTc37 = (10^(-7.54)) / (OneLiter/AmountOfSolutionIn1L);

      parameter Real L=L_old
        *
        (((KTh37/((10^(-7.2))+KTh37)) / (KRh37/((10^(-7.2))+KRh37)))^4)
        *
        (((KTz37*((10^(-7.2))^2 + KRz37*(10^(-7.2)) + KRz37*KRc37*(2.4217e-5)))/(KRz37*((10^(-7.2))^2 + KTz37*(10^(-7.2)) + KTz37*KTc37*(2.4217e-5))))^4)
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";

      Chemical.Components.Solution solution(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-56},{100,32}})));

      Chemical.Components.Reaction quaternaryForm(KC=KC)
        annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

      Chemical.Components.Substance CO_free(amountOfSubstance_start=initialCO*
            AmountOfSolutionIn1L, substanceData=Chemical.Examples.Substances.CarbonMonoxide_aqueous)
        annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
      Modelica.Blocks.Sources.Clock oxygenSource(offset=1)
        annotation (Placement(transformation(extent={{-78,48},{-58,68}})));
      Chemical.Sources.ExternalIdealGasSubstance CO_in_air(
        usePartialPressureInput=true,
        substanceData=Chemical.Examples.Substances.CarbonMonoxide_gas,
        Temperature=310.15) annotation (Placement(transformation(extent={{-10,-10},{
                10,10}}, origin={-36,58})));
      Chemical.Components.GasSolubility partialPressure1(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={-14,32})));

      Chemical.Components.Substance H2O(substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=AmountOfSolutionIn1L - THb - (initialO2 +
            initialCO2)*AmountOfSolutionIn1L)
        annotation (Placement(transformation(extent={{56,-54},{76,-34}})));
      HemoglobinQuaternaryForm_                             relaxed(
        DfG_selectedForm=DfG_tR,
        initialHb=THb/(L + 1),
        initialH=initialH,
        Kh37=KRh37,
        Kz37=KRz37,
        Kc37=KRc37,
        initialCO2=initialCO2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = -41000,
        Hz=8000,
        Hh=127000,
        initialCO=initialCO,
        Kco37=KRco37)
        annotation (Placement(transformation(extent={{-54,-44},{-34,-24}})));
      HemoglobinQuaternaryForm_                             tensed(
        DfG_selectedForm=DfG_tT,
        initialHb=THb*L/(L + 1),
        initialH=initialH,
        Kh37=KTh37,
        Kz37=KTz37,
        Kc37=KTc37,
        initialCO2=initialCO2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = 59000,
        Hz=-51000,
        Hh=59000,
        initialCO=initialCO,
        Kco37=KTco37)
                  annotation (Placement(transformation(extent={{32,-44},{12,-24}})));
      Chemical.Sources.ExternalMoleFraction H(substanceData=Chemical.Examples.Substances.Proton_aqueous,
          MoleFraction=initialH,
        Temperature=310.15)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}},
            origin={-12,-18})));
      Chemical.Components.Substance CO2_free(amountOfSubstance_start=initialCO2
            *AmountOfSolutionIn1L, substanceData(DfG_25degC_1bar=DfG_CO2,
            DfH_25degC=-412900))
        annotation (Placement(transformation(extent={{86,-8},{66,12}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
        PartialPressure(displayUnit="kPa") = pCO2,
        Temperature=310.15)                                   annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={62,60})));
      Chemical.Components.GasSolubility partialPressure2(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={62,32})));

    //  Real sO2 "Hemoglobin oxygen saturation";
      Real sCO "Hemoglobin carbon monoxide saturation";
      Real sCO2 "Hemoglobin carbon dioxide saturation";
      Real dH "Hemoglobin charge change caused by binding of Bohr's protons";
    equation
    //  sO2 = (sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x)) /
    //  (sum(relaxed.DeoxyHm.x) + sum(tensed.DeoxyHm.x) + sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x));

      sCO = (sum(relaxed.COHm.x) + sum(tensed.COHm.x)) /
      (sum(relaxed.DeoxyHm.x) + sum(tensed.DeoxyHm.x) + sum(relaxed.COHm.x) + sum(tensed.COHm.x));
    //+ sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x));

      sCO2 = (sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x)) /
      (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) + sum(relaxed.HmNH2.x) + sum(tensed.HmNH2.x) + sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x));

      dH = (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) - sum(relaxed.HmNHCOO.x) - sum(tensed.HmNHCOO.x) - sum(relaxed.HmA.x) - sum(tensed.HmA.x)) /
      THb;

      connect(oxygenSource.y, CO_in_air.partialPressure)
        annotation (Line(points={{-57,58},{-46,58}}, color={0,0,127}));

      connect(CO_in_air.port_a, partialPressure1.gas_port) annotation (Line(points={
              {-26,58},{-26,58},{-14,58},{-14,42}}, color={158,66,200}));
      connect(partialPressure1.liquid_port,CO_free. port_a) annotation (Line(
            points={{-14,22},{-14,-2},{-56,-2}}, color={158,66,200}));

      connect(CO_free.solution, solution.solution) annotation (Line(points={{-72,-12},
              {-72,-54},{60,-54},{60,-55.12}},        color={127,127,0}));
      connect(solution.solution, H2O.solution) annotation (Line(
          points={{60,-55.12},{60,-54}},
          color={127,127,0}));

      connect(relaxed.solution, solution.solution) annotation (Line(
          points={{-48,-42},{-48,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(relaxed.selectedForm, quaternaryForm.substrates[1]) annotation (Line(
          points={{-40,-42},{-22,-42}},
          color={158,66,200}));
      connect(tensed.solution, solution.solution) annotation (Line(
          points={{26,-42},{26,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(tensed.selectedForm, quaternaryForm.products[1]) annotation (Line(
          points={{18,-42},{-2,-42}},
          color={158,66,200}));
      connect(H.port_a, relaxed.H) annotation (Line(
          points={{-22,-18},{-36,-18},{-36,-26}},
          color={158,66,200}));
      connect(H.port_a, tensed.H) annotation (Line(
          points={{-22,-18},{14,-18},{14,-26}},
          color={158,66,200}));
      connect(CO2_gas.port_a, partialPressure2.gas_port) annotation (Line(
          points={{62,50},{62,42}},
          color={158,66,200}));
      connect(partialPressure2.liquid_port, CO2_free.port_a)
        annotation (Line(points={{62,22},{62,2},{66,2}}, color={158,66,200}));
      connect(CO2_free.port_a, tensed.CO2) annotation (Line(
          points={{66,2},{20,2},{20,-26}},
          color={158,66,200}));
      connect(CO2_free.port_a, relaxed.CO2) annotation (Line(
          points={{66,2},{-42,2},{-42,-26}},
          color={158,66,200}));
      connect(CO2_free.solution, solution.solution) annotation (Line(
          points={{82,-8},{82,-8},{82,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(CO_free.port_a, relaxed.CO) annotation (Line(
          points={{-56,-2},{-52,-2},{-52,-26}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(CO_free.port_a, tensed.CO) annotation (Line(
          points={{-56,-2},{30,-2},{30,-26}},
          color={158,66,200},
          smooth=Smooth.None));
      annotation (          experiment(StopTime=50),
        Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<p>Oxygen dissociation curve of hemoglobin.</p>
<p>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/ODC_pH_CO2_Hb.png\"/></p>
<p>J. W. Severinghaus, &quot;Simple, accurate equations for human blood O2 dissociation computations,&quot; Journal of Applied Physiology, vol. 46, pp. 599-602, 1979.</p>
<p><br>pO2 .. partial pressure of oxygen in gas</p>
<p>pCO2 .. partial pressure of carbon dioxide</p>
<p>sO2 .. oxygen saturation of hemoglobin</p>
<p>pH = log10(aH), where aH is mole fraction based activity of hydrogen ions</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/ODC_T_Hb.png\"/></p>
<p>R. B. Reeves, &quot;The effect of temperature on the oxygen equilibrium curve of human blood,&quot; Respiration physiology, vol. 42, pp. 317-328, 1980.</p>
<p><br>T .. temperature</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics),
        __Dymola_experimentSetupOutput);
    end HemoglobinMultipleAllostery_;

    model Joels57

      parameter Real r=5.42e-12, t=1.26e-9;
      HemoglobinMultipleAllostery_ pCO2_2000(pH=7.6, pCO2(displayUnit="kPa") = 2000,
        KRco37=r,
        KTco37=t)
        annotation (Placement(transformation(extent={{-50,-18},{-30,2}})));
      HemoglobinMultipleAllostery_ pCO2_5330(pH=7.5, pCO2(displayUnit="kPa") = 5330,
        KRco37=r,
        KTco37=t)
        annotation (Placement(transformation(extent={{-6,-16},{14,4}})));
      HemoglobinMultipleAllostery_ pCO2_9330(pH=7.4, pCO2(displayUnit="kPa") = 9330,
        KRco37=r,
        KTco37=t)
        annotation (Placement(transformation(extent={{40,-16},{60,4}})));
    end Joels57;

    model HemoglobinTitration "Multiple-ligand allosteric hemoglobin model"
      extends Modelica.Icons.Example;

      constant Modelica.SIunits.AmountOfSubstance THb = 0.001
        "Total amount of hemoglobin";

      constant Real RT=Modelica.Constants.R*298.15;

      constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        "Amount of solution used for molarity to mole fraction conversion";

      constant Modelica.SIunits.Volume OneLiter = 0.001;

      parameter Real L_old=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.SIunits.Concentration KR=0.000671946
        "oxygen dissociation on relaxed(R) hemoglobin subunit";

      parameter Modelica.SIunits.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.SIunits.MoleFraction KRo37 = KR*OneLiter/AmountOfSolutionIn1L;
      parameter Modelica.SIunits.MoleFraction KTo37 = KT*OneLiter/AmountOfSolutionIn1L;

      parameter Modelica.SIunits.MoleFraction KRco37 = (KR/218)*OneLiter/AmountOfSolutionIn1L;
      parameter Modelica.SIunits.MoleFraction KTco37 = (KT/218)*OneLiter/AmountOfSolutionIn1L;

      parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508);
      parameter Modelica.SIunits.ChemicalPotential DfG_CO2 = -RT*log(0.034/55.508)  + Chemical.Examples.Substances.CarbonDioxide_gas.DfG_25degC_1bar;

      parameter Modelica.SIunits.ChemicalPotential DfG_tT = 0;
      parameter Modelica.SIunits.ChemicalPotential DfG_tR = DfG_tT + RT * log(L);

      parameter Real KC = 1e-2 "Slow down factor";

      parameter Modelica.SIunits.MoleFraction initialO2=1.9594e-07
        "Initial O2 at 37degC, pO2=100Pa";             //at 25degC: 2.342e-8;
      parameter Modelica.SIunits.MoleFraction initialH=10^(-7.2);
      parameter Modelica.SIunits.MoleFraction initialCO2=2.4217e-10
        "Initial CO2 at 37degC, pCO2=40mmHg";          //at 25degC: 3.267e-5;

      parameter Modelica.SIunits.MoleFraction KRh37 = 10^(-6.89);
      parameter Modelica.SIunits.MoleFraction KTh37 = 10^(-7.52);

      parameter Modelica.SIunits.MoleFraction KRz37 = 10^(-7.25);
      parameter Modelica.SIunits.MoleFraction KTz37 = 10^(-7.73);

      parameter Modelica.SIunits.MoleFraction KRc37 = (10^(-8.35)) / (OneLiter/AmountOfSolutionIn1L);
      parameter Modelica.SIunits.MoleFraction KTc37 = (10^(-7.54)) / (OneLiter/AmountOfSolutionIn1L);

      parameter Real L=L_old
        *
        (((KTh37/((10^(-7.2))+KTh37)) / (KRh37/((10^(-7.2))+KRh37)))^4)
        *
        (((KTz37*((10^(-7.2))^2 + KRz37*(10^(-7.2)) + KRz37*KRc37*(2.4217e-5)))/(KRz37*((10^(-7.2))^2 + KTz37*(10^(-7.2)) + KTz37*KTc37*(2.4217e-5))))^4)
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";

      Chemical.Components.Solution solution(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-28},{100,60}})));

      Chemical.Components.Reaction quaternaryForm(KC=KC)
        annotation (Placement(transformation(extent={{-22,-24},{-2,-4}})));

      Chemical.Components.Substance O2_free(amountOfSubstance_start=initialO2*
            AmountOfSolutionIn1L, substanceData(DfG_25degC_1bar=DfG_O2,
            DfH_25degC=-11700))
        annotation (Placement(transformation(extent={{-76,28},{-56,48}})));
      Modelica.Blocks.Sources.Clock pHSource(offset=6.9)
        annotation (Placement(transformation(extent={{-60,-60},{-40,-40}})));
      Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
                                        substanceData=Chemical.Examples.Substances.Oxygen_gas,
        PartialPressure(displayUnit="bar") = 100000,
        Temperature=310.15)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            origin={-36,86})));
      Chemical.Components.GasSolubility partialPressure1(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={-14,60})));

      Chemical.Components.Substance H2O(substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=AmountOfSolutionIn1L - THb - (initialO2 +
            initialCO2)*AmountOfSolutionIn1L)
        annotation (Placement(transformation(extent={{56,-26},{76,-6}})));
      HemoglobinQuaternaryForm                              relaxed(
        Ko37=KRo37,
        DfG_selectedForm=DfG_tR,
        initialO2=initialO2,
        initialHb=THb/(L + 1),
        initialH=initialH,
        Kh37=KRh37,
        Kz37=KRz37,
        Kc37=KRc37,
        initialCO2=initialCO2,
        DfG_O2=DfG_O2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = -41000,
        Hz=8000,
        Hh=127000,
        Kco37=KRco37)
        annotation (Placement(transformation(extent={{-54,-16},{-34,4}})));
      HemoglobinQuaternaryForm                              tensed(
        Ko37=KTo37,
        DfG_selectedForm=DfG_tT,
        initialO2=initialO2,
        initialHb=THb*L/(L + 1),
        initialH=initialH,
        Kh37=KTh37,
        Kz37=KTz37,
        Kc37=KTc37,
        initialCO2=initialCO2,
        DfG_O2=DfG_O2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = 59000,
        Hz=-51000,
        Hh=59000,
        Kco37=KTco37)
                  annotation (Placement(transformation(extent={{32,-16},{12,4}})));
      Chemical.Sources.ExternalMoleFraction H(substanceData=Chemical.Examples.Substances.Proton_aqueous,
          MoleFraction=initialH,
        useMoleFractionInput=true,
        Temperature=310.15)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}},
            origin={-10,14})));
      Chemical.Components.Substance CO2_free(amountOfSubstance_start=initialCO2
            *AmountOfSolutionIn1L, substanceData(DfG_25degC_1bar=DfG_CO2,
            DfH_25degC=-412900))
        annotation (Placement(transformation(extent={{76,20},{56,40}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
        PartialPressure(displayUnit="kPa") = 5330,
        Temperature=310.15)                                   annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            origin={34,84})));
      Chemical.Components.GasSolubility partialPressure2(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={48,60})));

      Modelica.Blocks.Math.Exp exp
        annotation (Placement(transformation(extent={{4,-60},{24,-40}})));
      Modelica.Blocks.Math.Gain gain(k=-log(10))
        annotation (Placement(transformation(extent={{-28,-60},{-8,-40}})));

      Real sO2 "Hemoglobin oxygen saturation";
      Real sCO2 "Hemoglobin carbon dioxide saturation";
      Real dH "Hemoglobin charge change caused by binding of Bohr's protons";
    equation
      sO2 = (sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x)) /
      (sum(relaxed.DeoxyHm.x) + sum(tensed.DeoxyHm.x) + sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x));

      sCO2 = (sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x)) /
      (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) + sum(relaxed.HmNH2.x) + sum(tensed.HmNH2.x) + sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x));

      dH = (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) - sum(relaxed.HmNHCOO.x) - sum(tensed.HmNHCOO.x) - sum(relaxed.HmA.x) - sum(tensed.HmA.x)) /
      THb;

      connect(oxygen_in_air.port_a, partialPressure1.gas_port) annotation (Line(
          points={{-26,86},{-26,86},{-14,86},{-14,70}},
          color={158,66,200}));
      connect(partialPressure1.liquid_port, O2_free.port_a) annotation (Line(
            points={{-14,50},{-14,38},{-56,38}}, color={158,66,200}));

      connect(O2_free.solution, solution.solution) annotation (Line(points={{-72,28},
              {-72,-26},{60,-26},{60,-27.12}},        color={127,127,0}));
      connect(solution.solution, H2O.solution) annotation (Line(
          points={{60,-27.12},{60,-26}},
          color={127,127,0}));

      connect(relaxed.solution, solution.solution) annotation (Line(
          points={{-48,-14},{-48,-26},{60,-26},{60,-27.12}},
          color={127,127,0}));
      connect(relaxed.O2, O2_free.port_a) annotation (Line(
          points={{-52,2},{-52,38},{-56,38}},
          color={158,66,200}));
      connect(relaxed.selectedForm, quaternaryForm.substrates[1]) annotation (Line(
          points={{-40,-14},{-22,-14}},
          color={158,66,200}));
      connect(tensed.solution, solution.solution) annotation (Line(
          points={{26,-14},{26,-26},{60,-26},{60,-27.12}},
          color={127,127,0}));
      connect(tensed.O2, O2_free.port_a) annotation (Line(
          points={{30,2},{30,38},{-56,38}},
          color={158,66,200}));
      connect(tensed.selectedForm, quaternaryForm.products[1]) annotation (Line(
          points={{18,-14},{-2,-14}},
          color={158,66,200}));
      connect(H.port_a, relaxed.H) annotation (Line(
          points={{-20,14},{-36,14},{-36,2}},
          color={158,66,200}));
      connect(H.port_a, tensed.H) annotation (Line(
          points={{-20,14},{-24,14},{-24,2},{14,2}},
          color={158,66,200}));
      connect(CO2_gas.port_a, partialPressure2.gas_port) annotation (Line(
          points={{44,84},{44,84},{48,84},{48,70}},
          color={158,66,200}));
      connect(partialPressure2.liquid_port, CO2_free.port_a)
        annotation (Line(points={{48,50},{48,30},{56,30}},
                                                         color={158,66,200}));
      connect(CO2_free.port_a, tensed.CO2) annotation (Line(
          points={{56,30},{20,30},{20,2}},
          color={158,66,200}));
      connect(CO2_free.port_a, relaxed.CO2) annotation (Line(
          points={{56,30},{-42,30},{-42,2}},
          color={158,66,200}));
      connect(CO2_free.solution, solution.solution) annotation (Line(
          points={{72,20},{82,20},{82,-26},{60,-26},{60,-27.12}},
          color={127,127,0}));
      connect(exp.y, H.moleFractionInput) annotation (Line(
          points={{25,-50},{42,-50},{42,14},{0,14}},
          color={0,0,127}));
      connect(pHSource.y, gain.u) annotation (Line(
          points={{-39,-50},{-30,-50}},
          color={0,0,127}));
      connect(exp.u, gain.y) annotation (Line(
          points={{2,-50},{-7,-50}},
          color={0,0,127}));
      annotation (          experiment(StopTime=0.6),
        Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<p>Bohr&apos;s protons released during hemoglobin oxygenation</p>
<p>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
<p><br><img src=\"modelica://Chemical/Resources/Images/Examples/Bohr_Hb.png\"/></p>
<p>O. Siggaard-Andersen, &quot;Oxygen-Linked Hydrogen Ion Binding of Human Hemoglobin. Effects of Carbon Dioxide and 2, 3-Diphosphoglycerate I. Studies on Erythrolysate,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, vol. 27, pp. 351-360, 1971.</p>
</html>"));
    end HemoglobinTitration;

    model HemoglobinCarboxylation "Multiple-ligand allosteric hemoglobin model"
      extends Modelica.Icons.Example;

      constant Modelica.SIunits.AmountOfSubstance THb = 0.001
        "Total amount of hemoglobin";

      constant Real RT=Modelica.Constants.R*298.15;

      constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        "Amount of solution used for molarity to mole fraction conversion";

      constant Modelica.SIunits.Volume OneLiter = 0.001;

      parameter Real L_old=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.SIunits.Concentration KR=0.000671946
        "oxygen dissociation on relaxed(R) hemoglobin subunit";

      parameter Modelica.SIunits.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.SIunits.MoleFraction KRo37 = KR*OneLiter/AmountOfSolutionIn1L;
      parameter Modelica.SIunits.MoleFraction KTo37 = KT*OneLiter/AmountOfSolutionIn1L;

      parameter Modelica.SIunits.MoleFraction KRco37 = (KR/218)*OneLiter/AmountOfSolutionIn1L;
      parameter Modelica.SIunits.MoleFraction KTco37 = (KT/218)*OneLiter/AmountOfSolutionIn1L;

      parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508);
      parameter Modelica.SIunits.ChemicalPotential DfG_CO2 = -RT*log(0.034/55.508)  + Chemical.Examples.Substances.CarbonDioxide_gas.DfG_25degC_1bar;

      parameter Modelica.SIunits.ChemicalPotential DfG_tT = 0;
      parameter Modelica.SIunits.ChemicalPotential DfG_tR = DfG_tT + RT * log(L);

      parameter Real KC = 1e-3 "Slow down factor";

      parameter Modelica.SIunits.MoleFraction initialO2=1.9594e-07
        "Initial O2 at 37degC, pO2=100Pa";             //at 25degC: 2.342e-8;
      parameter Modelica.SIunits.MoleFraction initialH=10^(-7.2);
      parameter Modelica.SIunits.MoleFraction initialCO2=2.4217e-10
        "Initial CO2 at 37degC, pCO2=40mmHg";          //at 25degC: 3.267e-5;

      parameter Modelica.SIunits.MoleFraction KRh37 = 10^(-6.89);
      parameter Modelica.SIunits.MoleFraction KTh37 = 10^(-7.52);

      parameter Modelica.SIunits.MoleFraction KRz37 = 10^(-7.25);
      parameter Modelica.SIunits.MoleFraction KTz37 = 10^(-7.73);

      parameter Modelica.SIunits.MoleFraction KRc37 = (10^(-8.35)) / (OneLiter/AmountOfSolutionIn1L);
      parameter Modelica.SIunits.MoleFraction KTc37 = (10^(-7.54)) / (OneLiter/AmountOfSolutionIn1L);

      parameter Real L=L_old
        *
        (((KTh37/((10^(-7.2))+KTh37)) / (KRh37/((10^(-7.2))+KRh37)))^4)
        *
        (((KTz37*((10^(-7.2))^2 + KRz37*(10^(-7.2)) + KRz37*KRc37*(2.4217e-5)))/(KRz37*((10^(-7.2))^2 + KTz37*(10^(-7.2)) + KTz37*KTc37*(2.4217e-5))))^4)
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";

      Chemical.Components.Solution solution(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-56},{100,32}})));

      Chemical.Components.Reaction quaternaryForm(KC=KC)
        annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

      Chemical.Components.Substance O2_free(amountOfSubstance_start=initialO2*
            AmountOfSolutionIn1L, substanceData(DfG_25degC_1bar=DfG_O2,
            DfH_25degC=-11700))
        annotation (Placement(transformation(extent={{-76,0},{-56,20}})));
      Modelica.Blocks.Sources.Clock pHSource(offset=6.9)
        annotation (Placement(transformation(extent={{-60,-88},{-40,-68}})));
      Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
                                        substanceData=Chemical.Examples.Substances.Oxygen_gas,
        PartialPressure(displayUnit="bar") = 100000,
        Temperature=310.15)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            origin={-36,58})));
      Chemical.Components.GasSolubility partialPressure1(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={-14,32})));

      Chemical.Components.Substance H2O(substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=AmountOfSolutionIn1L - THb - (initialO2 +
            initialCO2)*AmountOfSolutionIn1L)
        annotation (Placement(transformation(extent={{56,-54},{76,-34}})));
      HemoglobinQuaternaryForm                              relaxed(
        Ko37=KRo37,
        DfG_selectedForm=DfG_tR,
        initialO2=initialO2,
        initialHb=THb/(L + 1),
        initialH=initialH,
        Kh37=KRh37,
        Kz37=KRz37,
        Kc37=KRc37,
        initialCO2=initialCO2,
        DfG_O2=DfG_O2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = -41000,
        Hz=8000,
        Hh=127000,
        Kco37=KRco37)
        annotation (Placement(transformation(extent={{-54,-44},{-34,-24}})));
      HemoglobinQuaternaryForm                              tensed(
        Ko37=KTo37,
        DfG_selectedForm=DfG_tT,
        initialO2=initialO2,
        initialHb=THb*L/(L + 1),
        initialH=initialH,
        Kh37=KTh37,
        Kz37=KTz37,
        Kc37=KTc37,
        initialCO2=initialCO2,
        DfG_O2=DfG_O2,
        DfG_CO2=DfG_CO2,
        KC=KC,
        Hc(displayUnit="kJ/mol") = 59000,
        Hz=-51000,
        Hh=59000,
        Kco37=KTco37)
                  annotation (Placement(transformation(extent={{32,-44},{12,-24}})));
      Chemical.Sources.ExternalMoleFraction H(substanceData=Chemical.Examples.Substances.Proton_aqueous,
          MoleFraction=initialH,
        useMoleFractionInput=true,
        Temperature=310.15)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}},
            origin={-10,-14})));
      Chemical.Components.Substance CO2_free(amountOfSubstance_start=initialCO2
            *AmountOfSolutionIn1L, substanceData(DfG_25degC_1bar=DfG_CO2,
            DfH_25degC=-412900))
        annotation (Placement(transformation(extent={{76,-8},{56,12}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
        PartialPressure(displayUnit="kPa") = 5330,
        Temperature=310.15)                                   annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            origin={34,56})));
      Chemical.Components.GasSolubility partialPressure2(useWaterCorrection=
            false, KC=KC) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={48,32})));

      Modelica.Blocks.Math.Exp exp
        annotation (Placement(transformation(extent={{4,-88},{24,-68}})));
      Modelica.Blocks.Math.Gain gain(k=-log(10))
        annotation (Placement(transformation(extent={{-28,-88},{-8,-68}})));

      Real sO2 "Hemoglobin oxygen saturation";
      Real sCO2 "Hemoglobin carbon dioxide saturation";
      Real dH "Hemoglobin charge change caused by binding of Bohr's protons";
    equation
      sO2 = (sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x)) /
      (sum(relaxed.DeoxyHm.x) + sum(tensed.DeoxyHm.x) + sum(relaxed.OxyHm.x) + sum(tensed.OxyHm.x));

      sCO2 = (sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x)) /
      (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) + sum(relaxed.HmNH2.x) + sum(tensed.HmNH2.x) + sum(relaxed.HmNHCOO.x) + sum(tensed.HmNHCOO.x));

      dH = (sum(relaxed.HmNH3.x) + sum(tensed.HmNH3.x) - sum(relaxed.HmNHCOO.x) - sum(tensed.HmNHCOO.x) - sum(relaxed.HmA.x) - sum(tensed.HmA.x)) /
      THb;

      connect(oxygen_in_air.port_a, partialPressure1.gas_port) annotation (Line(
          points={{-26,58},{-26,58},{-14,58},{-14,42}},
          color={158,66,200}));
      connect(partialPressure1.liquid_port, O2_free.port_a) annotation (Line(
            points={{-14,22},{-14,10},{-56,10}}, color={158,66,200}));

      connect(O2_free.solution, solution.solution) annotation (Line(points={{
              -72,0},{-72,-54},{60,-54},{60,-55.12}}, color={127,127,0}));
      connect(solution.solution, H2O.solution) annotation (Line(
          points={{60,-55.12},{60,-54}},
          color={127,127,0}));

      connect(relaxed.solution, solution.solution) annotation (Line(
          points={{-48,-42},{-48,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(relaxed.O2, O2_free.port_a) annotation (Line(
          points={{-52,-26},{-52,10},{-56,10}},
          color={158,66,200}));
      connect(relaxed.selectedForm, quaternaryForm.substrates[1]) annotation (Line(
          points={{-40,-42},{-22,-42}},
          color={158,66,200}));
      connect(tensed.solution, solution.solution) annotation (Line(
          points={{26,-42},{26,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(tensed.O2, O2_free.port_a) annotation (Line(
          points={{30,-26},{30,10},{-56,10}},
          color={158,66,200}));
      connect(tensed.selectedForm, quaternaryForm.products[1]) annotation (Line(
          points={{18,-42},{-2,-42}},
          color={158,66,200}));
      connect(H.port_a, relaxed.H) annotation (Line(
          points={{-20,-14},{-36,-14},{-36,-26}},
          color={158,66,200}));
      connect(H.port_a, tensed.H) annotation (Line(
          points={{-20,-14},{-24,-14},{-24,-26},{14,-26}},
          color={158,66,200}));
      connect(CO2_gas.port_a, partialPressure2.gas_port) annotation (Line(
          points={{44,56},{44,56},{48,56},{48,42}},
          color={158,66,200}));
      connect(partialPressure2.liquid_port, CO2_free.port_a)
        annotation (Line(points={{48,22},{48,2},{56,2}}, color={158,66,200}));
      connect(CO2_free.port_a, tensed.CO2) annotation (Line(
          points={{56,2},{20,2},{20,-26}},
          color={158,66,200}));
      connect(CO2_free.port_a, relaxed.CO2) annotation (Line(
          points={{56,2},{-42,2},{-42,-26}},
          color={158,66,200}));
      connect(CO2_free.solution, solution.solution) annotation (Line(
          points={{72,-8},{82,-8},{82,-54},{60,-54},{60,-55.12}},
          color={127,127,0}));
      connect(exp.y, H.moleFractionInput) annotation (Line(
          points={{25,-78},{42,-78},{42,-14},{0,-14}},
          color={0,0,127}));
      connect(pHSource.y, gain.u) annotation (Line(
          points={{-39,-78},{-30,-78}},
          color={0,0,127}));
      connect(exp.u, gain.y) annotation (Line(
          points={{2,-78},{-7,-78}},
          color={0,0,127}));
      annotation (          experiment(StopTime=0.6),
        Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<p>Binding of carbon dioxide to hemoglobin</p>
<p><br>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
<p><br><img src=\"modelica://Chemical/Resources/Images/Examples/CDDC_Hb.png\"/></p>
<p>C. Bauer and E. Schr&ouml;der, &quot;Carbamino compounds of haemoglobin in human adult and foetal blood,&quot; Journal of Physiology, vol. 227, pp. 457-471, 1972.</p>
</html>"));
    end HemoglobinCarboxylation;
  end Hemoglobin;

  package CheckSubstancesData
    model SimpleReaction
      "The simple chemical reaction A<->B with equilibrium B/A = 2"
       extends Modelica.Icons.Example;

      constant Real K = 2 "Dissociation constant of the reaction";

      constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

    Chemical.Sensors.DissociationCoefficient dissociationCoefficient
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    Chemical.Sources.PureSubstance A
        annotation (Placement(transformation(extent={{-56,-10},{-36,10}})));
    Chemical.Sources.PureSubstance B(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible, substanceData(DfG_25degC_1bar=-
              R*T_25degC*log(K)))
        annotation (Placement(transformation(extent={{60,-10},{40,10}})));
    equation
    connect(A.port_a, dissociationCoefficient.substrates[1]) annotation (Line(
        points={{-36,0},{-10,0}},
        color={158,66,200},
        thickness=1));
    connect(B.port_a, dissociationCoefficient.products[1]) annotation (Line(
        points={{40,0},{10,0}},
        color={158,66,200},
        thickness=1));
      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.001));
    end SimpleReaction;

    model SimpleReaction2
      "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
       extends Modelica.Icons.Example;

      constant Real Kb(unit="kg/mol") = 2
        "Molarity based dissociation constant of the reaction with one more reactant";

      constant Real Kx(unit="1") = Kb*55.508
        "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

      constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

    Chemical.Sources.PureSubstance A
        annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
      Chemical.Sensors.DissociationCoefficient reaction(nS=2)
        annotation (Placement(transformation(extent={{4,-8},{24,12}})));
    Chemical.Sources.PureSubstance B
        annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
    Chemical.Sources.PureSubstance C(substanceData(DfG_25degC_1bar=-R*T_25degC*
              log(Kx)))
        annotation (Placement(transformation(extent={{68,-8},{48,12}})));

    equation
      connect(reaction.products[1], C.port_a) annotation (Line(
          points={{24,2},{48,2}},
          color={158,66,200},
          thickness=1));

      connect(B.port_a, reaction.substrates[1]) annotation (Line(
          points={{-14,-14},{-10,-14},{-10,1.5},{4,1.5}},
          color={158,66,200},
          thickness=1));
      connect(A.port_a, reaction.substrates[2]) annotation (Line(
          points={{-14,12},{-10,12},{-10,2.5},{4,2.5}},
          color={158,66,200},
          thickness=1));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.001));
    end SimpleReaction2;

    model SimpleReaction2_Get_DfG
      "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
       extends Modelica.Icons.Example;

    Chemical.Sources.PureSubstance A
        annotation (Placement(transformation(extent={{-28,42},{-8,62}})));
      Chemical.Sensors.DissociationCoefficient reaction(nS=2)
        annotation (Placement(transformation(extent={{10,32},{30,52}})));
    Chemical.Sources.PureSubstance B
        annotation (Placement(transformation(extent={{-28,16},{-8,36}})));

      Modelica.Blocks.Math.InverseBlockConstraints inverseBlockConstraints
        annotation (Placement(transformation(extent={{-42,-80},{82,80}})));
    Chemical.Sources.ExternalElectroChemicalPotential C(usePotentialInput=true)
        annotation (Placement(transformation(extent={{60,32},{40,52}})));
      Modelica.Blocks.Sources.Constant K(k=2*55.508)
        annotation (Placement(transformation(extent={{-92,-10},{-72,10}})));
    equation

      connect(B.port_a, reaction.substrates[1]) annotation (Line(
          points={{-8,26},{-4,26},{-4,41.5},{10,41.5}},
          color={158,66,200},
          thickness=1));
      connect(A.port_a, reaction.substrates[2]) annotation (Line(
          points={{-8,52},{-4,52},{-4,42.5},{10,42.5}},
          color={158,66,200},
          thickness=1));
      connect(reaction.DissociationCoefficient_MoleFractionBased,
        inverseBlockConstraints.u2) annotation (Line(
          points={{20,34},{20,0},{-29.6,0}},
          color={0,0,127}));
      connect(reaction.products[1], C.port_a) annotation (Line(
          points={{30,42},{40,42}},
          color={158,66,200},
          thickness=1));
      connect(C.uInput, inverseBlockConstraints.y2) annotation (Line(
          points={{60,42},{70,42},{70,24},{46,24},{46,0},{72.7,0}},
          color={0,0,127}));
      connect(inverseBlockConstraints.u1, K.y) annotation (Line(
          points={{-48.2,0},{-71,0}},
          color={0,0,127}));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.001));
    end SimpleReaction2_Get_DfG;

    model StandardElectrochemicalCell
      "Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential "
     extends Modelica.Icons.Example;

      Chemical.Components.Solution cathode(ElectricGround=false)
        annotation (Placement(transformation(extent={{-90,-40},{-46,68}})));

      Chemical.Components.Solution anode(ElectricGround=false)
        annotation (Placement(transformation(extent={{60,-40},{96,70}})));

    Chemical.Sources.PureSubstance Ag(substanceData=Substances.Silver_solid)
        annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
    Chemical.Sources.PureSubstance Cl(substanceData=Substances.Chloride_aqueous)
        annotation (Placement(transformation(extent={{-8,-36},{-28,-16}})));
    Chemical.Sources.PureSubstance AgCl(substanceData=Substances.SilverChloride_solid)
        annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
    Chemical.Sources.ExternalIdealGasSubstance H2(
        substanceData=Substances.Hydrogen_gas,
        PartialPressure=100000,
        TotalPressure=100000)
        annotation (Placement(transformation(extent={{24,32},{44,52}})));
    Chemical.Sources.PureSubstance H(substanceData=Substances.Proton_aqueous)
        annotation (Placement(transformation(extent={{18,-36},{38,-16}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Chemical.Components.Reaction electrodeReaction(nP=2, p={2,2}) annotation (
         Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={52,6})));
      Chemical.Components.Reaction electrodeReaction1(nS=2, nP=2) annotation (
          Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-40,6})));
    Chemical.Components.ElectronTransfer electrone
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    Chemical.Components.ElectronTransfer electrone1
        annotation (Placement(transformation(extent={{86,-26},{66,-6}})));

    equation
      connect(Ag.port_a, electrodeReaction1.substrates[1]) annotation (Line(
          points={{-60,-18},{-42,-18},{-42,-4},{-42,-4}},
          color={158,66,200},
          thickness=1));
      connect(Cl.port_a, electrodeReaction1.substrates[2]) annotation (Line(
          points={{-28,-26},{-38,-26},{-38,-4}},
          color={158,66,200},
          thickness=1));
      connect(AgCl.port_a, electrodeReaction1.products[1]) annotation (Line(
          points={{-60,22},{-42,22},{-42,16}},
          color={158,66,200},
          thickness=1));
      connect(H2.port_a, electrodeReaction.substrates[1]) annotation (Line(
          points={{44,42},{52,42},{52,16}},
          color={158,66,200},
          thickness=1));
      connect(H.port_a, electrodeReaction.products[1]) annotation (Line(
          points={{38,-26},{54,-26},{54,-4}},
          color={158,66,200},
          thickness=1));
      connect(electrodeReaction.products[2], electrone1.port_a) annotation (Line(
          points={{50,-4},{50,-16},{66,-16}},
          color={158,66,200},
          thickness=1));
      connect(electrodeReaction1.products[2], electrone.port_a) annotation (Line(
          points={{-38,16},{-38,50},{-60,50}},
          color={158,66,200},
          thickness=1));
      connect(electrone.pin, voltageSensor.p) annotation (Line(
          points={{-80,50},{-86,50},{-86,74},{-6,74}},
          color={0,0,255}));
      connect(electrone1.pin, voltageSensor.n) annotation (Line(
          points={{86,-16},{90,-16},{90,74},{14,74}},
          color={0,0,255}));
    connect(electrone1.solution, anode.solution) annotation (Line(
        points={{82,-26},{80,-26},{80,-38.9},{88.8,-38.9}},
        color={127,127,0}));
    connect(electrone.solution, cathode.solution) annotation (Line(
        points={{-76,40},{-88,40},{-88,-38.92},{-54.8,-38.92}},
        color={127,127,0}));
      annotation (
      experiment(StopTime=1), Documentation(info=
                    "<html>
<p>Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential </p>
</html>", revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end StandardElectrochemicalCell;

    model StandardLeadAcidPotential
      "Standard potential of the lead acid battery"
     extends Modelica.Icons.Example;

      Chemical.Components.Solution anode(ElectricGround=false)
        annotation (Placement(transformation(extent={{54,-46},{92,62}})));

      Chemical.Components.Solution cathode(ElectricGround=false)
        annotation (Placement(transformation(extent={{-94,-50},{-56,58}})));

    Chemical.Sources.PureSubstance Pb(substanceData=Substances.Lead_solid)
        annotation (Placement(transformation(extent={{84,-34},{64,-14}})));
    Chemical.Sources.PureSubstance HSO4(substanceData=Substances.HydrogenSulfate_aqueous)
        annotation (Placement(transformation(extent={{-22,-58},{-2,-38}})));
    Chemical.Sources.PureSubstance PbSO4_(substanceData=Substances.LeadSulfate_solid)
        annotation (Placement(transformation(extent={{84,4},{64,24}})));
    Chemical.Sources.PureSubstance H(substanceData=Substances.Proton_aqueous)
        annotation (Placement(transformation(extent={{6,-28},{26,-8}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,60},{14,80}})));
      Chemical.Components.Reaction electrodeReaction(
        nP=2,
        nS=4,
        s={1,1,3,2},
        p={1,2}) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-42,14})));
      Chemical.Components.Reaction electrodeReaction1(
        nS=2,
        nP=3,
        p={1,1,2}) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={44,14})));

    Chemical.Components.ElectronTransfer electrone
        annotation (Placement(transformation(extent={{84,32},{64,52}})));
    Chemical.Components.ElectronTransfer electrone1
        annotation (Placement(transformation(extent={{-86,-12},{-66,8}})));
    Chemical.Sources.PureSubstance PbO2(substanceData=Substances.LeadDioxide_solid)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin=
               {-74,-30})));
    Chemical.Sources.PureSubstance H2O(substanceData=Substances.Water_liquid)
        annotation (Placement(transformation(extent={{-2,-10},{-22,10}})));
    Chemical.Sources.PureSubstance PbSO4(substanceData=Substances.LeadSulfate_solid)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin=
               {-74,32})));

    equation
      connect(Pb.port_a, electrodeReaction1.substrates[1]) annotation (Line(
          points={{64,-24},{45.5,-24},{45.5,4},{46,4}},
          color={158,66,200},
          thickness=1));
      connect(HSO4.port_a, electrodeReaction1.substrates[2]) annotation (Line(
          points={{-2,-48},{44,-48},{44,4},{42,4}},
          color={158,66,200},
          thickness=1));
      connect(PbSO4_.port_a, electrodeReaction1.products[1]) annotation (Line(
          points={{64,14},{56,14},{56,28},{46,28},{46,24},{46.6667,24}},
          color={158,66,200},
          thickness=1));
      connect(electrodeReaction.products[1], PbSO4.port_a) annotation (Line(
          points={{-44,24},{-44,32},{-64,32}},
          color={158,66,200},
          thickness=1));
      connect(electrodeReaction.products[2], H2O.port_a) annotation (Line(
          points={{-40,24},{-40,24},{-40,32},{-34,32},{-34,0},{-22,0}},
          color={158,66,200},
          thickness=1));
      connect(PbO2.port_a, electrodeReaction.substrates[1]) annotation (Line(
          points={{-64,-30},{-42,-30},{-42,4},{-45,4}},
          color={158,66,200},
          thickness=1));
      connect(HSO4.port_a, electrodeReaction.substrates[2]) annotation (Line(
          points={{-2,-48},{-40,-48},{-40,4},{-43,4}},
          color={158,66,200},
          thickness=1));
      connect(H.port_a, electrodeReaction.substrates[3]) annotation (Line(
          points={{26,-18},{-38,-18},{-38,4},{-41,4}},
          color={158,66,200},
          thickness=1));
      connect(electrone1.port_a, electrodeReaction.substrates[4]) annotation (Line(
          points={{-66,-2},{-44,-2},{-44,4},{-39,4}},
          color={158,66,200},
          thickness=1));
      connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(
          points={{26,-18},{32,-18},{32,32},{42,32},{42,24},{44,24}},
          color={158,66,200},
          thickness=1));
      connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
          points={{64,42},{44,42},{44,24},{41.3333,24}},
          color={158,66,200},
          thickness=1));
    connect(electrone1.pin, voltageSensor.p) annotation (Line(
        points={{-86,-2},{-98,-2},{-98,70},{-6,70}},
        color={0,0,255}));
    connect(electrone.pin, voltageSensor.n) annotation (Line(
        points={{84,42},{96,42},{96,70},{14,70}},
        color={0,0,255}));
    connect(electrone1.solution, cathode.solution) annotation (Line(
        points={{-82,-12},{-82,-48.92},{-63.6,-48.92}},
        color={127,127,0}));
    connect(electrone.solution, anode.solution) annotation (Line(
        points={{80,32},{80,-44.92},{84.4,-44.92}},
        color={127,127,0}));
      annotation (
      experiment(StopTime=100), Documentation(revisions=
                      "<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end StandardLeadAcidPotential;

  end CheckSubstancesData;

  model FluidAdapter
   extends Modelica.Icons.Example;

   package Medium = Chemical.Interfaces.SimpleChemicalMedium(substanceNames={"H2O(l)"},substanceData={Substances.Water_liquid});

    inner Modelica.Fluid.System system
      annotation (Placement(transformation(extent={{-82,66},{-62,86}})));
    Chemical.Components.FluidAdapter fluidConversion1(
      substanceNames={"H2O(l)"},
      substanceData={Substances.Water_liquid},
      redeclare package Medium = Medium)
      annotation (Placement(transformation(extent={{-50,-2},{-30,18}})));
    Chemical.Components.Solution simpleSolution1(BasePressure=110000,
        useThermalPort=true)
      annotation (Placement(transformation(extent={{-96,-20},{-26,40}})));
    Chemical.Components.Substance H2O1(
      redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
      substanceData=Substances.Water_liquid,
      amountOfSubstance_start=55)
      annotation (Placement(transformation(extent={{-80,-2},{-60,18}})));
    Chemical.Components.Solution simpleSolution2(temperature_start=299.15,
        useThermalPort=true)
      annotation (Placement(transformation(extent={{24,-20},{98,42}})));
    Chemical.Components.Substance H2O2(
      redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
      substanceData=Substances.Water_liquid,
      amountOfSubstance_start=55)
      annotation (Placement(transformation(extent={{84,-2},{64,18}})));
    Chemical.Components.FluidAdapter fluidConversion2(
      substanceNames={"H2O(l)"},
      substanceData={Substances.Water_liquid},
      redeclare package Medium = Medium)
      annotation (Placement(transformation(extent={{56,-2},{36,18}})));
    Modelica.Fluid.Pipes.StaticPipe pipe1(
      length=1,
      diameter=0.005,
      redeclare package Medium = Medium)
      annotation (Placement(transformation(extent={{-10,-2},{10,18}})));
  equation
  connect(fluidConversion1.solution, simpleSolution1.solution) annotation (
      Line(
      points={{-36,5},{-36,-8},{-40,-8},{-40,-19.4}},
      color={127,127,0}));
  connect(H2O1.port_a, fluidConversion1.substances[1]) annotation (Line(
      points={{-60,8},{-50,8}},
      color={158,66,200}));
  connect(H2O1.solution, simpleSolution1.solution) annotation (Line(
      points={{-76,-2},{-76,-8},{-40,-8},{-40,-19.4}},
      color={127,127,0}));
  connect(fluidConversion2.solution, simpleSolution2.solution) annotation (
      Line(
      points={{42,5},{42,-8},{80,-8},{80,-20},{83.2,-20},{83.2,-19.38}},
      color={127,127,0}));
  connect(H2O2.solution, simpleSolution2.solution) annotation (Line(
      points={{80,-2},{80,-19.38},{83.2,-19.38}},
      color={127,127,0}));
  connect(H2O2.port_a, fluidConversion2.substances[1]) annotation (Line(
      points={{64,8},{56,8}},
      color={158,66,200}));
  connect(fluidConversion1.fluid, pipe1.port_a) annotation (Line(
      points={{-30,8},{-10,8}},
      color={0,127,255}));
  connect(fluidConversion2.fluid, pipe1.port_b) annotation (Line(
      points={{36,8},{10,8}},
      color={0,127,255}));
    connect(simpleSolution1.heatPort, fluidConversion1.heatPort) annotation (Line(
        points={{-82,-20.6},{-82,-6},{-44,-6},{-44,5}},
        color={191,0,0}));
    connect(simpleSolution2.heatPort, fluidConversion2.heatPort) annotation (Line(
        points={{38.8,-20.62},{38.8,-6},{50,-6},{50,5}},
        color={191,0,0}));
    annotation (    experiment(StopTime=31));
  end FluidAdapter;

  model FluidAdapter2
   extends Modelica.Icons.Example;

   package Medium = Chemical.Interfaces.SimpleChemicalMedium(substanceNames={"H2O(l)","Ethanol"},substanceData={Substances.Water_liquid,
            Substances.Ethanol_liquid});

    inner Modelica.Fluid.System system
      annotation (Placement(transformation(extent={{-82,66},{-62,86}})));
    Chemical.Components.FluidAdapter fluidConversion1(
    substanceData={Substances.Water_liquid,Substances.Ethanol_liquid},
    substanceNames={"H2O(l)","Ethanol"},
    redeclare package Medium = Medium)
      annotation (Placement(transformation(extent={{-50,-2},{-30,18}})));

    Chemical.Components.Solution simpleSolution1(BasePressure=110000,
        useThermalPort=true)
      annotation (Placement(transformation(extent={{-96,-20},{-26,40}})));
    Chemical.Components.Substance H2O(
      redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
      substanceData=Substances.Water_liquid,
      amountOfSubstance_start=5)
      annotation (Placement(transformation(extent={{-90,-2},{-70,18}})));
    Chemical.Components.Solution simpleSolution2(temperature_start=299.15,
        useThermalPort=true)
      annotation (Placement(transformation(extent={{24,-20},{98,42}})));
    Chemical.Components.Substance H2O_(
      redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
      substanceData=Substances.Water_liquid,
      amountOfSubstance_start=55)
      annotation (Placement(transformation(extent={{80,-2},{60,18}})));
    Chemical.Components.FluidAdapter fluidConversion2(
    substanceNames={"H2O(l)","Ethanol"},
    substanceData={Substances.Water_liquid,Substances.Ethanol_liquid},
    redeclare package Medium =Medium)
      annotation (Placement(transformation(extent={{56,-2},{36,18}})));

    Modelica.Fluid.Pipes.StaticPipe pipe1(
      length=1,
    diameter=0.005,
    redeclare package Medium =Medium)
                  annotation (Placement(transformation(extent={{-10,-2},{10,
            18}})));
    Chemical.Components.Substance C2H5OH(
    redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
    substanceData=Substances.Ethanol_liquid,
      amountOfSubstance_start=100)
    annotation (Placement(transformation(extent={{-70,18},{-50,38}})));
    Chemical.Components.Substance C2H5OH_(
    redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
    amountOfSubstance_start=55,
    substanceData=Substances.Ethanol_liquid)
    annotation (Placement(transformation(extent={{90,20},{70,40}})));
  equation
  connect(fluidConversion1.solution, simpleSolution1.solution) annotation (
      Line(
      points={{-36,5},{-36,-8},{-40,-8},{-40,-19.4}},
      color={127,127,0}));
    connect(H2O.port_a, fluidConversion1.substances[1]) annotation (Line(
        points={{-70,8},{-60,8},{-60,6},{-50,6}},
        color={158,66,200}));
    connect(H2O.solution, simpleSolution1.solution) annotation (Line(
        points={{-86,-2},{-86,-8},{-40,-8},{-40,-19.4}},
        color={127,127,0}));
  connect(fluidConversion2.solution, simpleSolution2.solution) annotation (
      Line(
      points={{42,5},{42,-6},{82,-6},{82,-19.38},{83.2,-19.38}},
      color={127,127,0}));
  connect(H2O_.solution, simpleSolution2.solution) annotation (Line(
      points={{76,-2},{76,-6},{82,-6},{82,-19.38},{83.2,-19.38}},
      color={127,127,0}));
  connect(H2O_.port_a, fluidConversion2.substances[1]) annotation (Line(
      points={{60,8},{58,8},{58,6},{56,6}},
      color={158,66,200}));
  connect(fluidConversion1.fluid, pipe1.port_a) annotation (Line(
      points={{-30,8},{-10,8}},
      color={0,127,255}));
  connect(fluidConversion2.fluid, pipe1.port_b) annotation (Line(
      points={{36,8},{10,8}},
      color={0,127,255}));
    connect(simpleSolution1.heatPort, fluidConversion1.heatPort) annotation (Line(
        points={{-82,-20.6},{-82,-6},{-44,-6},{-44,5}},
        color={191,0,0}));
    connect(simpleSolution2.heatPort, fluidConversion2.heatPort) annotation (Line(
        points={{38.8,-20.62},{38.8,-4},{50,-4},{50,5}},
        color={191,0,0}));
  connect(C2H5OH.solution, simpleSolution1.solution) annotation (Line(
      points={{-66,18},{-66,-8},{-40,-8},{-40,-19.4}},
      color={127,127,0}));
  connect(C2H5OH.port_a, fluidConversion1.substances[2]) annotation (Line(
      points={{-50,28},{-50,10}},
      color={158,66,200}));
  connect(C2H5OH_.solution, simpleSolution2.solution) annotation (Line(
      points={{86,20},{86,-6},{82,-6},{82,-19.38},{83.2,-19.38}},
      color={127,127,0}));
  connect(C2H5OH_.port_a, fluidConversion2.substances[2]) annotation (Line(
      points={{70,30},{56,30},{56,10}},
      color={158,66,200}));
    annotation (    experiment(
        StopTime=18.4),
      Documentation(info="<html>
<p>Demonstration of compatibility with FluidPort from Modelica Standard Library.</p>
</html>"));
  end FluidAdapter2;
end Examples;
