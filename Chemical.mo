within ;
package Chemical
  "Chemical library (reactions, diffusions, semipermeable membranes, gas dissolutions, electrochemical cells, ...)"
 extends Modelica.Icons.Package;

  package Examples
    "Examples that demonstrate usage of the Pressure flow components"
  extends Modelica.Icons.ExamplesPackage;

    package Substances "Definitions of substances"
        extends Modelica.Icons.Package;

      constant Interfaces.SubstanceModel.SubstanceData Silver_solid(
        MolarWeight=0.1078682,
        z=0,
        DfH=0,
        DfG_25degC=0,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"}) "Ag(s)";

      constant Interfaces.SubstanceModel.SubstanceData Silver_aqueous(
        MolarWeight=0.1078682,
        z=1,
        DfH=105900,
        DfG_25degC=77100,
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "Ag+(aq)";

      constant Interfaces.SubstanceModel.SubstanceData SilverChloride_solid(
        MolarWeight=0.14332,
        z=0,
        DfH=-127030,
        DfG_25degC=-109720,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "AgCl(s)";

      constant Interfaces.SubstanceModel.SubstanceData Calcium_aqueous(
        MolarWeight=0.0401,
        z=2,
        DfH=-542960,
        DfG_25degC=-542960 - 298.15*(33.67),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
        "}) "Ca++(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Chloride_aqueous(
        MolarWeight=0.03545,
        z=-1,
        DfH=-167460,
        DfG_25degC=-131170,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "Cl-(aq)";

      constant Interfaces.SubstanceModel.SubstanceData CarbonDioxide_gas(
        MolarWeight=0.044,
        DfH=-393500,
        DfG_25degC=-394400,
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "CO2(g)";

      constant Interfaces.SubstanceModel.SubstanceData CarbonDioxide_aqueous(
        MolarWeight=0.044,
        DfH=-412900,
        DfG_25degC=-386200,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "CO2(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Carbonate_aqueous(
        MolarWeight=0.06001,
        z=-2,
        DfH=-676300,
        DfG_25degC=-676300 - 298.15*(-497.065),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "CO3--(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Electrone_solid(
        MolarWeight=5.4857990946e-7,
        z=-1,
        DfH=0,
        DfG_25degC=0,
        References={"http://physics.nist.gov/cgi-bin/cuu/Value?mme",
            "To solve standard electo-chemical cell potentials"}) "e-(s)";

      constant Interfaces.SubstanceModel.SubstanceData Iron2_aqueous(
        MolarWeight=0.05585,
        z=2,
        DfH=-87860,
        DfG_25degC=-87860 - 298.15*(-9.93),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
        "}) "Fe++(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Iron3_aqueous(
        MolarWeight=0.05585,
        z=3,
        DfH=-47700,
        DfG_25degC=-47700 - 298.15*(-124.77),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
"}) "Fe+++(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Glucose_solid(
        MolarWeight=0.1806,
        DfH=-1274500,
        DfG_25degC=-1274500 - 298.15*(-1220.66),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
"}) "Glu(s)";

      constant Interfaces.SubstanceModel.SubstanceData Hydrogen_gas(
        MolarWeight=0.00201588,
        z=0,
        DfH=0,
        DfG_25degC=0,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"}) "H2(g)";

      constant Interfaces.SubstanceModel.SubstanceData CarbonicAcid_aqueous(
        MolarWeight=0.062027,
        DfH=-699700,
        DfG_25degC=-699700 - 298.15*(-256.582),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H2CO3(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Water_gas(
        MolarWeight=0.018015,
        DfH=-241830,
        DfG_25degC=-228590,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H2O(g)";

      constant Interfaces.SubstanceModel.SubstanceData Water_liquid(
        MolarWeight=0.018015,
        DfH=-285830,
        DfG_25degC=-237190,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H2O(l)";

      constant Interfaces.SubstanceModel.SubstanceData
        DihydrogenPhosphate_aqueous(
        MolarWeight=0.095,
        z=-1,
        DfH=-1302480,
        DfG_25degC=-1302480 - 298.15*(-561.395),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "H2PO4-(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Hydronium_aqueous(
        MolarWeight=0.019022,
        z=1,
        DfH=-285840,
        DfG_25degC=-285840 - 298.15*(-163.17),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H3O+(aq)";

      constant Interfaces.SubstanceModel.SubstanceData PhosphoricAcid_aqueous(
        MolarWeight=0.095,
        DfH=-1288000,
        DfG_25degC=-1288000 - 298.15*(-496.4),
        References={"https://en.wikipedia.org/wiki/Phosphoric_acid",
            "https://www.researchgate.net/publication/6600409_Standard_thermodynamic_properties_of_H3PO4%28aq%29_over_a_wide_range_of_temperatures_and_pressures"})
        "H3PO4(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Proton_aqueous(
        MolarWeight=0.001007,
        z=1,
        DfH=0,
        DfG_25degC=0,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H+(aq)";
                 // as hypothetical HA <-> H+ + A- simplification of H2O + HA <-> H3O+ + A-";

      constant Interfaces.SubstanceModel.SubstanceData Bicarbonate_aqueous(
        MolarWeight=0.06102,
        z=-1,
        DfH=-691100,
        DfG_25degC=-691100 - 298.15*(-348.82),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "HCO3-(aq)";

      constant Interfaces.SubstanceModel.SubstanceData
        HydrogenPhosphate_aqueous(
        MolarWeight=0.095,
        z=-2,
        DfH=-1298700,
        DfG_25degC=-1298700 - 298.15*(-686.232),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "HPO4--(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Potassium_aqueous(
        MolarWeight=0.0391,
        z=1,
        DfH=-251200,
        DfG_25degC=-251200 - 298.15*(103.97),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "K+(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Magnesium_aqueous(
        MolarWeight=0.0243,
        z=2,
        DfH=-461960,
        DfG_25degC=-461960 - 298.15*(-19.99),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
","http://www.vias.org/genchem/standard_enthalpies_table.html"}) "Mg++(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Sodium_aqueous(
        MolarWeight=0.02299,
        z=1,
        DfH=-239660,
        DfG_25degC=-239660 - 298.15*(74.49),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "Na+(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Amonium_aqueous(
        MolarWeight=0.01804,
        z=1,
        DfH=-132800,
        DfG_25degC=-132800 - 298.15*(-178.77),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
"}) "NH4+(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Oxygen_gas(
        MolarWeight=0.032,
        DfH=0,
        DfG_25degC=0,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"}) "O2(g)";

      constant Interfaces.SubstanceModel.SubstanceData Oxygen_aqueous(
        MolarWeight=0.032,
        DfH=-11700,
        DfG_25degC=16320,
        References={
            "http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.pdf",
            "https://books.google.cz/books?id=dr-VBAAAQBAJ&pg=PA156&lpg=PA156&dq=Gibbs+energy+of+formation++%22O2(aq)%22&source=bl&ots=09N5CxY7OD&sig=hbsTXQvX59vXBqHUjFVVIZQpHCA&hl=cs&sa=X&ei=sDQtVaeUMMaRsAHpzYHgAg&redir_esc=y#v=onepage&q=Gibbs%20energy%20of%20formation%20%20%22O2(aq)%22&f=false"})
        "O2(aq)";
      constant Interfaces.SubstanceModel.SubstanceData Hydroxide_aqueous(
        MolarWeight=0.017006,
        z=-1,
        DfH=-229940,
        DfG_25degC=-157300,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "OH-(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Phosphate_aqueous(
        MolarWeight=0.095,
        z=-3,
        DfH=-1284070,
        DfG_25degC=-1284070 - 298.15*(-866.946),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "PO4---(aq)";

      constant Interfaces.SubstanceModel.SubstanceData Sulphates_aqueous(
        MolarWeight=0.09607,
        z=-2,
        DfH=-907500,
        DfG_25degC=-907500 - 298.15*(-555.123),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "SO4--(aq)";

        //Some organic molecules: https://www.e-education.psu.edu/drupal6/files/be497b/pdf/Bioenergetics_AppA.pdf
    end Substances;

    model SimpleReaction
      "The simple chemical reaction A<->B with equilibrium B/A = 2"
       extends Modelica.Icons.Example;

      constant Real K = 2 "Dissociation constant of the reaction";

      constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Components.Substance A(
        amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-64,-8},{-44,12}})));

      Components.Reaction reaction annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
      Components.Substance B(
        substanceData( DfG_25degC=-R*T_25degC*log(K)),
        amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{42,-8},{62,12}})));
      Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    equation

      connect(reaction.products[1], B.port_a) annotation (Line(
          points={{10,2},{62,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(A.solution, solution.solution) annotation (Line(
          points={{-60,-8},{-60,-92},{0,-92},{0,-100}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(B.solution, solution.solution) annotation (Line(points={{46,-8},{
              46,-92},{0,-92},{0,-100}}, smooth=Smooth.None));
      connect(A.port_a, reaction.substrates[1]) annotation (Line(
          points={{-44,2},{-10,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.001),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics),
        __Dymola_experimentSetupOutput);
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

      Components.Substance A(amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
      Components.Reaction reaction(nS=2)
        annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Components.Substance B(amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Components.Substance C(substanceData(DfG_25degC=-R*T_25degC*log(Kx)), amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{68,-8},{48,12}})));
      Components.Solution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    equation

      connect(reaction.products[1], C.port_a) annotation (Line(
          points={{24,2},{24,2},{48,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(A.solution, solution.solution) annotation (Line(
          points={{-30,2},{-30,-90},{0,-90},{0,-100}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(C.solution, solution.solution) annotation (Line(points={{64,-8},{66,-8},
              {66,-90},{0,-90},{0,-100}},        smooth=Smooth.None));
      connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},{-30,
              -90},{0,-90},{0,-100}},    smooth=Smooth.None));

      connect(B.port_a, reaction.substrates[1]) annotation (Line(
          points={{-14,-14},{-10,-14},{-10,1.5},{4,1.5}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(A.port_a, reaction.substrates[2]) annotation (Line(
          points={{-14,12},{-10,12},{-10,2.5},{4,2.5}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.001),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics),
        __Dymola_experimentSetupOutput);
    end SimpleReaction2;

    model ExothermicReaction
      "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"

       extends Modelica.Icons.Example;

      parameter Modelica.SIunits.MolarEnergy ReactionEnthalpy=-55000;

      Components.Substance A( amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-56,-60},{-36,-40}})));
      Components.Reaction reaction
        annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
      Components.Substance B( amountOfSubstance_start=0.1, substanceData(DfH=ReactionEnthalpy))
        annotation (Placement(transformation(extent={{44,-60},{64,-40}})));
      Components.Solution thermal_isolated_solution(ConstantTemperature=false)
        annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));

      Modelica.SIunits.HeatFlowRate q
        "Heat flow to environment to reach constant temperature";
      Modelica.SIunits.Temperature t
        "Temperature if the solution is ideally thermal isolated from environment";

      Components.Substance A1(amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-56,40},{-36,60}})));
      Components.Reaction reaction1
        annotation (Placement(transformation(extent={{-8,40},{12,60}})));
      Components.Substance B1(amountOfSubstance_start=0.1, substanceData(DfH=-ReactionEnthalpy))
        annotation (Placement(transformation(extent={{44,40},{64,60}})));
      Components.Solution solution_at_constant_temperature
        annotation (Placement(transformation(extent={{-100,0},{98,94}})));
    equation
      q = -solution_at_constant_temperature.heatFromEnvironment;

      t = thermal_isolated_solution.solution.T;

      connect(A.port_a, reaction.substrates[1]) annotation (Line(
          points={{-36,-50},{-8,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(reaction.products[1], B.port_a) annotation (Line(
          points={{12,-50},{64,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(B.solution, thermal_isolated_solution.solution) annotation (Line(
          points={{48,-60},{48,-88},{-1,-88},{-1,-100}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(A.solution, thermal_isolated_solution.solution) annotation (Line(
            points={{-52,-60},{-52,-88},{-1,-88},{-1,-100}}, smooth=Smooth.None));
      connect(A1.port_a, reaction1.substrates[1]) annotation (Line(
          points={{-36,50},{-8,50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(reaction1.products[1], B1.port_a) annotation (Line(
          points={{12,50},{64,50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(B1.solution, solution_at_constant_temperature.solution) annotation (
          Line(
          points={{48,40},{48,12},{-1,12},{-1,0}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(A1.solution, solution_at_constant_temperature.solution) annotation (
          Line(points={{-52,40},{-52,12},{-1,12},{-1,0}}, smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.001),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                        graphics),
        __Dymola_experimentSetupOutput);
    end ExothermicReaction;

    model Henry "Dissolution of gases in liquids"
       extends Modelica.Icons.Example;
      Components.GasSolubility CO2_new
        annotation (Placement(transformation(extent={{-86,-20},{-66,0}})));
      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,

      Components.Substance CO2_l_plasma(substanceData=
            Substances.CarbonDioxide_aqueous)
        annotation (Placement(transformation(extent={{-86,-60},{-66,-40}})));
      Components.GasSolubility O2_new
        annotation (Placement(transformation(extent={{8,-20},{28,0}})));

      Sources.AirSubstance O2_g_n1(
        substanceData=Substances.Oxygen_gas,
        PartialPressure=12665.626804425,
        TotalPressure=101325.0144354)
        annotation (Placement(transformation(extent={{32,34},{52,54}})));
      Components.Substance O2_l_plasma(substanceData=
            Substances.Oxygen_aqueous)
        annotation (Placement(transformation(extent={{10,-60},{30,-40}})));
      Components.GasSolubility CO2_new1
        annotation (Placement(transformation(extent={{-50,-20},{-30,0}})));

      Sources.AirSubstance CO2_g_n2(
        substanceData=Substances.CarbonDioxide_gas,
        PartialPressure=5332.8954966,
        TotalPressure=101325.0144354)
        annotation (Placement(transformation(extent={{-80,34},{-60,54}})));

      Components.Substance CO2_l_erythrocyte(substanceData=
            Substances.CarbonDioxide_aqueous)
        annotation (Placement(transformation(extent={{-50,-62},{-30,-42}})));
      Components.GasSolubility O2_new1
        annotation (Placement(transformation(extent={{40,-18},{60,2}})));

      Components.Substance O2_l_erythrocyte(substanceData=
            Substances.Oxygen_aqueous)
        annotation (Placement(transformation(extent={{40,-60},{60,-40}})));
      Components.Solution plasma(amountOfSolution_start=52.3)
        annotation (Placement(transformation(extent={{-104,-100},{96,100}})));
      Components.Solution redCells(amountOfSolution_start=39.7)
        annotation (Placement(transformation(extent={{-92,-90},{108,110}})));
      Components.GasSolubility O2_new2(useWaterCorrection=false)
        annotation (Placement(transformation(extent={{74,-20},{94,0}})));
      Components.Substance O2_l_erythrocyte1(substanceData(DfG_25degC=-Modelica.Constants.R*298.15*log(0.0013*0.018)))
        annotation (Placement(transformation(extent={{74,-60},{94,-40}})));
    equation

      connect(CO2_g_n2.port_a, CO2_new.gas_port) annotation (Line(
          points={{-60,44},{-52,44},{-52,12},{-76,12},{-76,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(CO2_g_n2.port_a, CO2_new1.gas_port) annotation (Line(
          points={{-60,44},{-52,44},{-52,12},{-40,12},{-40,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(CO2_new.liquid_port, CO2_l_plasma.port_a) annotation (Line(
          points={{-76,-20},{-76,-50},{-66,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(CO2_new1.liquid_port, CO2_l_erythrocyte.port_a) annotation (Line(
          points={{-40,-20},{-40,-36},{-40,-52},{-30,-52}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(O2_g_n1.port_a, O2_new.gas_port) annotation (Line(
          points={{52,44},{60,44},{60,18},{18,18},{18,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(O2_g_n1.port_a, O2_new1.gas_port) annotation (Line(
          points={{52,44},{60,44},{60,18},{50,18},{50,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(O2_new.liquid_port, O2_l_plasma.port_a) annotation (Line(
          points={{18,-20},{18,-50},{30,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(O2_new1.liquid_port, O2_l_erythrocyte.port_a) annotation (Line(
          points={{50,-18},{50,-34},{50,-50},{60,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(CO2_l_plasma.solution, plasma.solution) annotation (Line(
          points={{-82,-60},{-82,-100},{-4,-100}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(O2_l_plasma.solution, plasma.solution) annotation (Line(
          points={{14,-60},{14,-100},{-4,-100}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(CO2_l_erythrocyte.solution, redCells.solution) annotation (Line(
          points={{-46,-62},{-46,-90},{8,-90}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(O2_l_erythrocyte.solution, redCells.solution) annotation (Line(
          points={{44,-60},{44,-90},{8,-90}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(O2_g_n1.port_a, O2_new2.gas_port) annotation (Line(
          points={{52,44},{60,44},{60,18},{84,18},{84,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(O2_new2.liquid_port, O2_l_erythrocyte1.port_a) annotation (Line(
          points={{84,-20},{84,-50},{94,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(O2_l_erythrocyte1.solution, redCells.solution) annotation (Line(
          points={{78,-60},{78,-90},{8,-90}},
          color={0,0,0},
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}),
                graphics),
        experiment(StopTime=1),
        Documentation(info="<html>
<p>Please note, that the total content of CO2 and O2 in blood plasma and erythrocytes must be determined by including bicarbonate and hemoglobin connected amounts.  </p>
</html>"));
    end Henry;

    model MichaelisMenten "Basic enzyme kinetics"
      extends Modelica.Icons.Example;

      //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
      Sources.AmbientMoleFraction
                              P(             substanceData(DfG_25degC=-Modelica.Constants.R*298.15*50),
          MoleFraction=0.1)
        annotation (Placement(transformation(extent={{92,-12},{72,8}})));
      Sources.AmbientMoleFraction
                              S(MoleFraction=Km)
        annotation (Placement(transformation(extent={{-88,-34},{-68,-14}})));

         parameter Modelica.SIunits.AmountOfSubstance tE=0.01
        "Total amount of enzyme";
         parameter Real k_cat(unit="1/s", displayUnit="1/min")= 1
        "Forward rate of second reaction";
         constant Modelica.SIunits.Concentration Km=0.1
        "Michaelis constant = substrate concentration at rate of half Vmax";

        parameter Modelica.SIunits.MolarFlowRate Vmax=tE*k_cat
        "Maximal molar flow";
        parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution= 55.508
        "Amount of solution used in kinetics";

          Components.Substance ES(substanceData(DfG_25degC=-Modelica.Constants.R*298.15*log(2/Km)),
          amountOfSubstance_start=tE/2)
            annotation (Placement(transformation(extent={{-12,-10},{8,10}})));
          Components.Substance E(amountOfSubstance_start=tE/2)
            annotation (Placement(transformation(extent={{-10,38},{10,58}})));
      Components.Reaction chemicalReaction(nS=2, ActivationEnergy=2*AmountOfSolution*Modelica.Constants.R
            *298.15*log(2)/Vmax,
        AmountOfSolution=AmountOfSolution)
        annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));

      Components.Reaction chemicalReaction1(nP=2, ActivationEnergy=2*AmountOfSolution*Modelica.Constants.R
            *298.15*(50 - log(2))/Vmax,
        AmountOfSolution=AmountOfSolution)
        annotation (Placement(transformation(extent={{24,-10},{44,10}})));

      Components.Solution solution annotation (Placement(transformation(extent={{-100,
                -100},{100,100}})));
    equation

         //Michaelis-Menton: v=((E.q_out.conc + ES.q_out.conc)*k_cat)*S.concentration/(Km+S.concentration);

      connect(S.port_a, chemicalReaction.substrates[1]) annotation (Line(
          points={{-68,-24},{-56,-24},{-56,-0.5},{-42,-0.5}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(chemicalReaction.products[1], ES.port_a) annotation (Line(
          points={{-22,0},{8,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(ES.port_a, chemicalReaction1.substrates[1]) annotation (Line(
          points={{8,0},{24,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(E.port_a, chemicalReaction.substrates[2]) annotation (Line(
          points={{10,48},{-52,48},{-52,0.5},{-42,0.5}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(E.port_a, chemicalReaction1.products[2]) annotation (Line(
          points={{10,48},{54,48},{54,0.5},{44,0.5}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(chemicalReaction1.products[1], P.port_a) annotation (Line(
          points={{44,-0.5},{58,-0.5},{58,-2},{72,-2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(E.solution, solution.solution) annotation (Line(
          points={{-6,38},{-18,38},{-18,-92},{0,-92},{0,-100}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(ES.solution, solution.solution)
        annotation (Line(points={{-8,-10},{-8,-92},{0,-92},{0,-100}}, smooth=Smooth.None));
          annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>To recalculate the enzyme kinetics from Michaelis-Menton parameters Km, tE a k_cat is selected the same half-rate of the reaction defined as:</p>
<p>E = ES = tE/2 .. the amount of free enzyme is the same as the amount of enzyme-substrate complexes</p>
<p>S = Km .. the amount of substrate is Km</p>
<p>r = Vmax/2 = tE*k_cat / 2 .. the rate of reaction is the half of maximal rate</p>
<p><br>Conversions of molar concentration to mole fraction (MM is molar mass of the solvent in solution -&GT; 55.508 kg/mol for water):</p>
<p>x(Km) = Km/MM</p>
<p>x(tE) = tE/MM</p>
<p>xS = S/MM = Km/MM</p>
<p><br>The new kinetics of the system defined as:</p>
<p>uS&deg; = DfG(S) = 0</p>
<p>uE&deg; = DfG(E) = 0</p>
<p>uES&deg; = DfG(ES) = DfG(S) + DfG(E) -R*T*ln(2/x(Km))</p>
<p>uP&deg; = DfG(P) </p>
<p><br>r = Vmax/2</p>
<p>r = -C1 * (uES&deg; - uE&deg; - uS&deg; + R*T*ln(xES/(xE*xS) ) = -C1 * (-R*T*ln(2/x(Km)) - R*T*ln(xS) ) = C1 * R * T * ln(2) </p>
<p>r = -C2 * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) = -C2 * (DfG(P) - uES&deg; + R*T*ln(xP) ) = C2 * (-DfG(P) - R * T * ln(2))</p>
<p><br>ActivationEnergy1 = AmountOfSolution/(Tau*(Vmax/2)) * R * T * ln(2) </p>
<p>ActivationEnergy2 = AmountOfSolution/(Tau*(Vmax/2)) * ( -DfG(P) - R * T * ln(2) ) </p>
<p><br>where</p>
<p>AmountOfSolution = MM = 55.508 (for water)</p>
<p>Tau = 1 s (just to be physical unit correct)</p>
<p>DfG(P) = -R*T*50 is Gibbs energy of formation of product (setting negative enough makes second reaction almost irreversible)</p>
</html>"),
        experiment,
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                        graphics),
        __Dymola_experimentSetupOutput);
    end MichaelisMenten;

    model StandardElectrochemicalCell
      "Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential "
     extends Modelica.Icons.Example;

      Sources.PureSubstance Ag(substanceData=
            Chemical.Examples.Substances.Silver_solid)
        annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
      Sources.PureSubstance Cl(substanceData=
            Chemical.Examples.Substances.Chloride_aqueous)
                                       annotation (Placement(transformation(extent={{-8,-36},{-28,-16}})));
      Sources.PureSubstance AgCl(substanceData=
            Chemical.Examples.Substances.SilverChloride_solid)
        annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
      Sources.AirSubstance H2(substanceData = Chemical.Examples.Substances.Hydrogen_gas,
        PartialPressure=100000,
        TotalPressure=100000) annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Sources.PureSubstance H(substanceData=
            Chemical.Examples.Substances.Proton_aqueous)
                                      annotation (Placement(transformation(extent={{18,-36},{38,-16}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Components.Reaction electrodeReaction(nP=2, p={2,2}) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={52,6})));
      Components.Reaction electrodeReaction1(nS=2, nP=2) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-40,6})));
      Sources.PureElectricParticle
                           electrone(substanceData=
            Chemical.Examples.Substances.Electrone_solid)
                                     annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
      Sources.PureElectricParticle
                           electrone1(substanceData=
            Chemical.Examples.Substances.Electrone_solid)
                                      annotation (Placement(transformation(extent={{86,-26},{66,-6}})));
    equation
      connect(Ag.port_a, electrodeReaction1.substrates[1]) annotation (Line(
          points={{-60,-18},{-42,-18},{-42,-4},{-40.5,-4}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Cl.port_a, electrodeReaction1.substrates[2]) annotation (Line(
          points={{-28,-26},{-39.5,-26},{-39.5,-4}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(AgCl.port_a, electrodeReaction1.products[1]) annotation (Line(
          points={{-60,22},{-40.5,22},{-40.5,16}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(H2.port_a, electrodeReaction.substrates[1]) annotation (Line(
          points={{44,42},{52,42},{52,16}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(H.port_a, electrodeReaction.products[1]) annotation (Line(
          points={{38,-26},{52.5,-26},{52.5,-4}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(electrodeReaction.products[2], electrone1.port_a) annotation (Line(
          points={{51.5,-4},{51.5,-16},{66,-16}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(electrodeReaction1.products[2], electrone.port_a) annotation (Line(
          points={{-39.5,16},{-39.5,50},{-60,50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(electrone.pin, voltageSensor.p) annotation (Line(
          points={{-80,50},{-88,50},{-88,74},{-6,74}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(electrone1.pin, voltageSensor.n) annotation (Line(
          points={{86,-16},{94,-16},{94,74},{14,74}},
          color={0,0,255},
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics), Documentation(info="<html>
<p>Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential </p>
</html>"));
    end StandardElectrochemicalCell;

    model ElectrochemicalCell
      "The electrochemical cell: Pt(s) | H2(g) | H+(aq), Cl−(aq) | AgCl(s) | Ag(s)"
     extends Modelica.Icons.Example;

      Sources.PureSubstance Ag(substanceData=
            Chemical.Examples.Substances.Silver_solid)
        annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
      Components.Substance Cl(substanceData=
            Chemical.Examples.Substances.Chloride_aqueous,
          amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-28,-36},{-8,-16}})));
      Sources.PureSubstance AgCl(substanceData=
            Chemical.Examples.Substances.SilverChloride_solid)
        annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
      Sources.AirSubstance H2(substanceData = Chemical.Examples.Substances.Hydrogen_gas,
        PartialPressure=100000,
        TotalPressure=100000) annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Components.Substance H(substanceData=
            Chemical.Examples.Substances.Proton_aqueous,
        amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{10,-36},{30,-16}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Components.Reaction electrodeReaction(nP=2, p={2,2}) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={52,6})));
      Components.Reaction electrodeReaction1(nS=2, nP=2) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-40,6})));
      Components.Solution solution1 annotation (Placement(transformation(extent={{-32,-46},{42,-6}})));
      Sources.PureElectricParticle
                           electrone(substanceData=
            Chemical.Examples.Substances.Electrone_solid)
                                     annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
      Sources.PureElectricParticle
                           electrone1(substanceData=
            Chemical.Examples.Substances.Electrone_solid)
                                      annotation (Placement(transformation(extent={{86,-26},{66,-6}})));
    equation
      connect(Ag.port_a, electrodeReaction1.substrates[1]) annotation (Line(
          points={{-60,-18},{-42,-18},{-42,-4},{-40.5,-4}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Cl.port_a, electrodeReaction1.substrates[2]) annotation (Line(
          points={{-8,-26},{-39.5,-26},{-39.5,-4}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(AgCl.port_a, electrodeReaction1.products[1]) annotation (Line(
          points={{-60,22},{-40.5,22},{-40.5,16}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(H2.port_a, electrodeReaction.substrates[1]) annotation (Line(
          points={{44,42},{52,42},{52,16}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(H.port_a, electrodeReaction.products[1]) annotation (Line(
          points={{30,-26},{52.5,-26},{52.5,-4}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(electrodeReaction.products[2], electrone1.port_a) annotation (Line(
          points={{51.5,-4},{51.5,-16},{66,-16}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(electrodeReaction1.products[2], electrone.port_a) annotation (Line(
          points={{-39.5,16},{-39.5,50},{-60,50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(electrone.pin, voltageSensor.p) annotation (Line(
          points={{-80,50},{-88,50},{-88,74},{-6,74}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(electrone1.pin, voltageSensor.n) annotation (Line(
          points={{86,-16},{94,-16},{94,74},{14,74}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(Cl.solution, solution1.solution) annotation (Line(
          points={{-24,-36},{-24,-40},{5,-40},{5,-46}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(H.solution, solution1.solution) annotation (Line(points={{14,-36},
              {14,-40},{5,-40},{5,-46}}, smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),
            graphics));
    end ElectrochemicalCell;

    package AcidBase

      model WaterSelfIonization "H2O  <->  OH-   +   H+ "
        import Chemical;
          extends Modelica.Icons.Example;
        Components.Substance H3O(
          amountOfSubstance_start=1e-7,   substanceData=
              Chemical.Examples.Substances.Hydronium_aqueous)
                                        annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={36,72})));
        Components.Substance OH(
          amountOfSubstance_start=1e-7,   substanceData=
              Chemical.Examples.Substances.Hydroxide_aqueous)
                                        annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={36,26})));
        Components.Substance H2O(
          amountOfSubstance_start=1,   substanceData=
              Chemical.Examples.Substances.Water_liquid)
                                     annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={-40,46})));
        Chemical.Components.Reaction waterDissociation(nP=2, s={2})
          annotation (Placement(transformation(extent={{-12,36},{8,56}})));
              Real pH, pH_;
        Components.Substance H_(
          amountOfSubstance_start=1e-7,   substanceData=
              Chemical.Examples.Substances.Proton_aqueous)       annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={34,-30})));
        Components.Substance OH_(
          amountOfSubstance_start=1e-7,   substanceData=
              Chemical.Examples.Substances.Hydroxide_aqueous)
                                        annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={34,-76})));
        Components.Substance H2O_(
          amountOfSubstance_start=1,   substanceData=
              Chemical.Examples.Substances.Water_liquid)
                                     annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={-40,-56})));
        Chemical.Components.Reaction waterDissociation_(nP=2)
          annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));
        Chemical.Components.Solution solution(amountOfSolution_start=1)
          annotation (Placement(transformation(extent={{-74,4},{74,98}})));
        Chemical.Components.Solution solution1(amountOfSolution_start=1)
          annotation (Placement(transformation(extent={{-76,-94},{72,0}})));
      equation
        pH = -log10(                H3O.x);
                    /*H3O.gamma * */
        pH_ = -log10(               H_.x);
                     /*H_.gamma * */

        connect(OH.port_a, waterDissociation.products[1]) annotation (Line(
            points={{46,26},{22,26},{22,45.5},{8,45.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(waterDissociation.products[2], H3O.port_a) annotation (Line(
            points={{8,46.5},{22,46.5},{22,72},{46,72}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.port_a, waterDissociation.substrates[1]) annotation (Line(
            points={{-30,46},{-12,46}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(OH_.port_a,waterDissociation_. products[1]) annotation (Line(
            points={{44,-76},{20,-76},{20,-56.5},{6,-56.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(waterDissociation_.products[2], H_.port_a) annotation (Line(
            points={{6,-55.5},{20,-55.5},{20,-30},{44,-30}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O_.port_a,waterDissociation_. substrates[1]) annotation (Line(
            points={{-30,-56},{-14,-56}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.solution, solution.solution) annotation (Line(
            points={{-46,36},{0,36},{0,4}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(OH.solution, solution.solution) annotation (Line(
            points={{30,16},{30,4},{0,4}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H3O.solution, solution.solution) annotation (Line(
            points={{30,62},{52,62},{52,4},{0,4}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H2O_.solution, solution1.solution) annotation (Line(
            points={{-46,-66},{-2,-66},{-2,-94}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(OH_.solution, solution1.solution) annotation (Line(
            points={{28,-86},{28,-94},{-2,-94}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H_.solution, solution1.solution) annotation (Line(
            points={{28,-40},{52,-40},{52,-94},{-2,-94}},
            color={0,0,0},
            smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>Self-ionization of water.</p>
<p>Ions difference (SID) in water causes the acidity/basicity, where pH = -log10(aH+). An activity of hydrogen ions aH+ is approximated with concentration (mol/l) of the oxonium cations H3O+.</p>
<pre><b>plotExpression(apply(-log10(WaterSelfIonization.H3O.solute)),&nbsp;false,&nbsp;&QUOT;pH&QUOT;,&nbsp;1);</b></pre>
<p><br>The titration slope der(pH)/der(SID)=1.48e+6 1/(mol/L) at pH=7.4.</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(__Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}),            graphics),
          __Dymola_experimentSetupOutput);
      end WaterSelfIonization;

      model CarbonDioxideInWater "CO2 as alone acid-base buffer"
        import Chemical;
          extends Modelica.Icons.Example;
        Components.Substance HCO3(  substanceData=
              Chemical.Examples.Substances.Bicarbonate_aqueous)
          annotation (Placement(transformation(extent={{-14,24},{6,44}})));
        Chemical.Components.Reaction HendersonHasselbalch(nP=2, nS=2,
          Tau=1200) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-58,22},{-38,42}})));
        Sources.AirSubstance CO2_gas(
            substanceData=
              Chemical.Examples.Substances.CarbonDioxide_gas,
          PartialPressure=5332.8954966,
          TotalPressure=101325.0144354) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-80,86})));
        Components.Substance H(  substanceData=
              Chemical.Examples.Substances.Proton_aqueous,
            amountOfSubstance_start=3e-8) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}}, origin={-4,-10})));
        Components.GasSolubility gasSolubility(Tau=1200)
          annotation (Placement(transformation(extent={{-90,52},{-70,72}})));
                                              /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
        Components.Substance CO2_liquid(  substanceData=
              Chemical.Examples.Substances.CarbonDioxide_aqueous)
          annotation (Placement(transformation(extent={{-90,22},{-70,42}})));
        Components.Substance CO3(  substanceData=
              Chemical.Examples.Substances.Carbonate_aqueous)
          annotation (Placement(transformation(extent={{68,32},{88,52}})));
        Chemical.Components.Reaction c2(nP=2, nS=1)
          "K=10^(-10.33 + 3), dH=14.9kJ/mol"
          annotation (Placement(transformation(extent={{20,24},{40,44}})));
        Chemical.Components.Substance H2O(  substanceData=
              Chemical.Examples.Substances.Water_liquid,
            amountOfSubstance_start=55.507)
          annotation (Placement(transformation(extent={{-90,-22},{-70,-2}})));
        Real pH;
        Chemical.Components.Solution solution
          annotation (Placement(transformation(extent={{-100,-100},{100,80}})));
      equation
        pH = -log10(                                H.x);
                    /*H.port_a.activityCoefficient*/

        connect(CO2_gas.port_a, gasSolubility.gas_port) annotation (Line(
            points={{-80,76},{-80,72}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(gasSolubility.liquid_port, CO2_liquid.port_a) annotation (Line(
            points={{-80,52},{-80,42},{-80,32},{-70,32}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch.products[1], H.port_a) annotation (Line(
            points={{-38,31.5},{-20,31.5},{-20,-10},{6,-10}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch.products[2], HCO3.port_a) annotation (Line(
            points={{-38,32.5},{-20,32.5},{-20,34},{6,34}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HCO3.port_a, c2.substrates[1]) annotation (Line(
            points={{6,34},{20,34}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(c2.products[1], H.port_a) annotation (Line(
            points={{40,33.5},{52,33.5},{52,-10},{6,-10}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(c2.products[2], CO3.port_a) annotation (Line(
            points={{40,34.5},{56,34.5},{56,42},{88,42}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.port_a, HendersonHasselbalch.substrates[2]) annotation (
            Line(
            points={{-70,32},{-70,32},{-70,32.5},{-58,32.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.port_a, HendersonHasselbalch.substrates[1]) annotation (Line(
            points={{-70,-12},{-64,-12},{-64,31.5},{-58,31.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.solution, solution.solution) annotation (Line(
            points={{-86,22},{-86,-72},{0,-72},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H2O.solution, solution.solution) annotation (Line(points={{-86,-22},{-86,
                -72},{0,-72},{0,-100}}, smooth=Smooth.None));
        connect(HCO3.solution, solution.solution) annotation (Line(points={{-10,24},{-10,
                -72},{0,-72},{0,-100}}, smooth=Smooth.None));
        connect(H.solution, solution.solution) annotation (Line(points={{-10,-20},{-10,
                -72},{0,-72},{0,-100}}, smooth=Smooth.None));
        connect(CO3.solution, solution.solution) annotation (Line(points={{72,32},{72,
                -72},{0,-72},{0,-100}}, smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>CO2 solution in water without any other acid-base buffers.</p>
<pre><b>plotExpression(apply(-log10(CarbonDioxideInWater.H3O.solute)),&nbsp;false,&nbsp;&QUOT;pH&QUOT;,&nbsp;1);</b></pre>
<p><br>Please note, that OH- (and CO3^-2) can be neglected from electroneutrality calculation, because of very small concentrations (in physiological pH) anyway. </p>
<p>And if SID&GT;0 then also H3O+ can be also neglected from electroneutrality, because only bicarbonate anions HCO3- (or CO3^-2) are needed there to balance the electroneutrality.</p>
<p><br>The partial pressure of CO2 in gas are input parameter. Outputs are an amount of free disolved CO2 in liquid and an amount of HCO3-.</p>
<p><br>The titration slope der(pH)/der(SID)=17.5 1/(mol/L) at pH=7.4 and pCO2=40 mmHg.</p>
<p><br>Molar heat of formation (aqueous):</p>
<p>CO2:        -413.5 kJ/mol  (gas: -393.5 kJ/mol )</p>
<p>H2O:        -285.8 kJ/mol</p>
<p>HCO3-:        -692.0 kJ/mol</p>
<p>CO3^-2:        -677.1 kJ/mol</p>
<p><br>Enthalphy of reaction H2O + CO2 &LT;-&GT; HCO3- + H+  :         7.3 kJ/mol</p>
<p>Enthalphy of reaction HCO3- &LT;-&GT; CO3^-2 + H+  :        14.9 kJ/mol</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.02, __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}),            graphics),
          __Dymola_experimentSetupOutput);
      end CarbonDioxideInWater;

      model Phosphate
        import Chemical;
          extends Modelica.Icons.Example;

        parameter Physiolibrary.Types.Concentration totalPO4=0.00115
          "Total phosphate concentration";

        Components.Substance H(
          amountOfSubstance_start=55.6*10^(-7.4),
            substanceData=Chemical.Examples.Substances.Proton_aqueous)
          "hydrogen ions activity" annotation (Placement(transformation(extent=
                  {{-10,-10},{10,10}}, origin={36,-12})));

        Components.Substance H3PO4(
          amountOfSubstance_start=1e-08,
            substanceData=
              Chemical.Examples.Substances.PhosphoricAcid_aqueous)
          annotation (Placement(transformation(extent={{-98,-58},{-78,-38}})));
        Components.Substance H2PO4(
          amountOfSubstance_start=0.0005,
            substanceData=
              Chemical.Examples.Substances.DihydrogenPhosphate_aqueous)
          annotation (Placement(transformation(extent={{-44,-58},{-24,-38}})));
        Components.Substance HPO4(
            substanceData=
              Chemical.Examples.Substances.HydrogenPhosphate_aqueous,
          amountOfSubstance_start=0.0006)
          annotation (Placement(transformation(extent={{16,-58},{36,-38}})));
        Components.Substance PO4(
            substanceData=Chemical.Examples.Substances.Phosphate_aqueous,
          amountOfSubstance_start=1e-08)
          annotation (Placement(transformation(extent={{72,-58},{92,-38}})));

        Chemical.Components.Reaction chemicalReaction(nP=2) "10^(-1.915 + 3)"
          annotation (Placement(transformation(extent={{-70,-58},{-50,-38}})));
        Chemical.Components.Reaction chemicalReaction1(nP=2) "10^(-6.66 + 3)"
          annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
        Chemical.Components.Reaction chemicalReaction2(nP=2) "10^(-11.78 + 3)"
          annotation (Placement(transformation(extent={{44,-58},{64,-38}})));

        Chemical.Components.Solution solution
          annotation (Placement(transformation(extent={{-100,-100},{100,26}})));
      equation
        connect(H3PO4.port_a, chemicalReaction.substrates[1]) annotation (Line(
            points={{-78,-48},{-70,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction.products[1], H2PO4.port_a) annotation (Line(
            points={{-50,-48.5},{-42,-48.5},{-42,-48},{-24,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2PO4.port_a, chemicalReaction1.substrates[1]) annotation (Line(
            points={{-24,-48},{-14,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[1], HPO4.port_a) annotation (Line(
            points={{6,-48.5},{16,-48.5},{16,-48},{36,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(HPO4.port_a, chemicalReaction2.substrates[1]) annotation (Line(
            points={{36,-48},{44,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction2.products[1], PO4.port_a) annotation (Line(
            points={{64,-48.5},{74,-48.5},{74,-48},{92,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction.products[2], H.port_a) annotation (Line(
            points={{-50,-47.5},{-44,-47.5},{-44,-32},{46,-32},{46,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[2], H.port_a) annotation (Line(
            points={{6,-47.5},{14,-47.5},{14,-32},{46,-32},{46,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction2.products[2], H.port_a) annotation (Line(
            points={{64,-47.5},{72,-47.5},{72,-32},{46,-32},{46,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H3PO4.solution, solution.solution) annotation (Line(
            points={{-94,-58},{-94,-88},{0,-88},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H2PO4.solution, solution.solution) annotation (Line(points={{-40,-58},
                {-40,-88},{0,-88},{0,-100}}, smooth=Smooth.None));
        connect(HPO4.solution, solution.solution) annotation (Line(points={{20,-58},{22,
                -58},{22,-88},{0,-88},{0,-100}}, smooth=Smooth.None));
        connect(PO4.solution, solution.solution) annotation (Line(points={{76,-58},{76,
                -88},{0,-88},{0,-100}}, smooth=Smooth.None));
        connect(H.solution, solution.solution) annotation (Line(points={{30,-22},{30,-88},
                {0,-88},{0,-100}}, smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>Henderson-Hasselbalch equation in ideal buffered solution, where pH remains constant.</p>
<p>The partial pressure of CO2 in gas are input parameter. Outputs are an amount of free disolved CO2 in liquid and an amount of HCO3-.</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.05),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics));
      end Phosphate;

      model AlbuminTitration "Figge-Fencl model (22. Dec. 2007)"
        extends Modelica.Icons.Example;

        Sources.AmbientMoleFraction   H(MoleFraction=1e-7)
          "hydrogen ions activity"                                        annotation (Placement(
              transformation(extent={{10,-10},{-10,10}}, origin={14,22})));

        constant Integer n=218 "Number of weak acid group in albumin molecule";
        constant Real pKAs[n]=cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
          "acid dissociation constants";
        constant Real K[n]=fill(10.0, n) .^ (-pKAs);
        constant Real DfG[n]= Modelica.Constants.R*(298.15)*log(K);

        Chemical.Components.Substance A[n](
          each amountOfSubstance_start=0.00033) "deprotonated acid groups"
          annotation (Placement(transformation(extent={{4,-16},{24,4}})));
        Chemical.Components.Reaction react[n](each nP=2)
          annotation (Placement(transformation(extent={{-44,-2},{-24,18}})));

        Chemical.Components.Substance HA[n](substanceData(DfG_25degC=DfG), each amountOfSubstance_start=
             0.00033) "protonated acid groups"
          annotation (Placement(transformation(extent={{-76,-2},{-56,18}})));

        Components.Solution solution
          annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
      equation
        connect(react.products[1], A.port_a) annotation (Line(
            points={{-24,7.5},{-12,7.5},{-12,-6},{24,-6}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        for i in 1:n loop
          connect(react[i].products[2], H.port_a) annotation (Line(
              points={{-24,8.5},{-14,8.5},{-14,22},{4,22}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(HA[i].solution, solution.solution) annotation (Line(
            points={{-72,-2},{-72,-86},{0,-86},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
          connect(A[i].solution, solution.solution) annotation (Line(
            points={{8,-16},{8,-86},{0,-86},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        end for;
        connect(HA.port_a, react.substrates[1]) annotation (Line(
            points={{-56,8},{-44,8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        annotation ( Documentation(revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",       info="<html>
<pre><b>plotExpression(apply(-log10(AlbuminTitration.H3O.solute)),&nbsp;false,&nbsp;&QUOT;pH&QUOT;,&nbsp;1);</b></pre>
<p>The titration slope der(pH)/der(SID)=185 1/(mol/L) at pH=7.4 and tAlb=0.66 mmol/l.</p>
<p><br>Data and model is described in</p>
<p><font style=\"color: #222222; \">Jame Figge: Role of non-volatile weak acids (albumin, phosphate and citrate). In: Stewart&apos;s Textbook of Acid-Base, 2nd Edition, John A. Kellum, Paul WG Elbers editors, &nbsp;AcidBase org, 2009, pp. 216-232.</font></p>
</html>"),experiment(
            StopTime=1e-005,
            __Dymola_fixedstepsize=5e-005,
            __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics),
          __Dymola_experimentSetupOutput);
      end AlbuminTitration;



      model CarbonDioxideInBlood
        import Chemical;
          extends Modelica.Icons.Example;
        Components.Substance HCO3(
          amountOfSubstance_start=0.024,
          substanceData=Chemical.Examples.Substances.Bicarbonate_aqueous)
          annotation (Placement(transformation(extent={{10,-10},{-10,10}},
              rotation=0,
              origin={34,32})));
        Chemical.Components.Reaction HendersonHasselbalch(nP=2, nS=2)
          "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-20,22},{0,42}})));
        Sources.AirSubstance CO2_gas(
          substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
          PartialPressure(displayUnit="mmHg") = 5332.8954966,
          TotalPressure(displayUnit="mmHg") = 101325.0144354)
                                        annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-80,82})));
        Chemical.Sources.AmbientMoleFraction  H(
          MoleFraction = 10^(-7.4),
          substanceData=Chemical.Examples.Substances.Proton_aqueous)
                                                        annotation (Placement(
              transformation(extent={{-10,-10},{10,10}}, origin={56,8},
              rotation=180)));
        Components.GasSolubility gasSolubility
          annotation (Placement(transformation(extent={{-90,46},{-70,66}})));
                                              /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
        Components.Substance CO2_liquid(
          amountOfSubstance_start=0.00123,
          substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous)
          annotation (Placement(transformation(extent={{-92,22},{-72,42}})));
        Components.Substance H2O(
          amountOfSubstance_start=52,
          substanceData=Chemical.Examples.Substances.Water_liquid)
          annotation (Placement(transformation(extent={{-64,6},{-44,26}})));
        Components.Substance HCO3_E(
          amountOfSubstance_start=0.0116,
          substanceData=Chemical.Examples.Substances.Bicarbonate_aqueous)
          annotation (Placement(transformation(extent={{46,-70},{26,-50}})));
        Chemical.Components.Reaction HendersonHasselbalch1(nP=2, nS=2)
          "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-16,-70},{4,-50}})));
        Components.Substance CO2_liquid_E(
          amountOfSubstance_start=0.00093,
          substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous)
          annotation (Placement(transformation(extent={{-94,-86},{-74,-66}})));
        Components.Substance H2O_E(
          amountOfSubstance_start=39.5,
          substanceData=Chemical.Examples.Substances.Water_liquid)
          annotation (Placement(transformation(extent={{-60,-62},{-40,-42}})));
        Components.Substance Cl_E(
          amountOfSubstance_start=0.0499,
          substanceData=Chemical.Examples.Substances.Chloride_aqueous)
          annotation (Placement(transformation(extent={{98,-94},{78,-74}})));
        Components.Substance Cl_P(
          amountOfSubstance_start=0.103,
          substanceData=Chemical.Examples.Substances.Chloride_aqueous)
          annotation (Placement(transformation(extent={{96,2},{76,22}})));

        Real pH_e,pH_p;
        Chemical.Components.Solution blood_erythrocytes(amountOfSolution_start=
              39.7) annotation (Placement(transformation(extent={{-100,-100},{
                  100,-32}})));
        Chemical.Components.Solution blood_plasma(amountOfSolution_start=52.3)
          annotation (Placement(transformation(extent={{-100,-16},{100,64}})));
        Components.Membrane membrane1
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-34,-26})));
        Components.Membrane membrane2
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={18,-28})));
        Components.Membrane membrane3
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={72,-28})));
        Components.Membrane membrane annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-68,-28})));
        Chemical.Sources.AmbientMoleFraction H_E(substanceData=Chemical.Examples.Substances.Proton_aqueous,
            MoleFraction=10^(-7.2)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              origin={58,-82},
              rotation=180)));
      equation
        pH_p = -log10(H.x);
        pH_e = -log10(H_E.x);
        connect(HendersonHasselbalch.products[1], HCO3.port_a) annotation (Line(
            points={{0,31.5},{6,31.5},{6,32},{24,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H.port_a, HendersonHasselbalch.products[2]) annotation (Line(
            points={{46,8},{12,8},{12,32.5},{0,32.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.port_a, HendersonHasselbalch.substrates[1]) annotation (
           Line(
            points={{-72,32},{-70,32},{-70,31.5},{-20,31.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.port_a, HendersonHasselbalch.substrates[2]) annotation (Line(
            points={{-44,16},{-34,16},{-34,32.5},{-20,32.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch1.products[1], HCO3_E.port_a) annotation (Line(
            points={{4,-60.5},{16,-60.5},{16,-60},{26,-60}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid_E.port_a, HendersonHasselbalch1.substrates[1]) annotation (
            Line(
            points={{-74,-76},{-30,-76},{-30,-60.5},{-16,-60.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O_E.port_a, HendersonHasselbalch1.substrates[2]) annotation (Line(
            points={{-40,-52},{-34,-52},{-34,-59.5},{-16,-59.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.solution, blood_plasma.solution) annotation (Line(
            points={{-88,22},{-90,22},{-90,-16},{0,-16}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H2O.solution, blood_plasma.solution) annotation (Line(points={{
                -60,6},{-60,-16},{0,-16}}, smooth=Smooth.None));
        connect(Cl_P.solution, blood_plasma.solution) annotation (Line(
            points={{92,2},{92,-16},{0,-16}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(CO2_liquid_E.solution, blood_erythrocytes.solution) annotation
          (Line(
            points={{-90,-86},{-90,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H2O_E.solution, blood_erythrocytes.solution) annotation (Line(
              points={{-56,-62},{-56,-100},{0,-100}}, smooth=Smooth.None));
        connect(Cl_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{94,-94},{94,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(HCO3_E.solution, blood_erythrocytes.solution) annotation (Line(
              points={{42,-70},{42,-100},{0,-100}}, smooth=Smooth.None));
        connect(gasSolubility.liquid_port, CO2_liquid.port_a) annotation (Line(
            points={{-80,46},{-80,32},{-72,32}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(membrane.port_b, CO2_liquid_E.port_a) annotation (Line(
            points={{-68,-38},{-68,-76},{-74,-76}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(membrane1.port_b, H2O_E.port_a) annotation (Line(
            points={{-34,-36},{-34,-52},{-40,-52}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(membrane1.port_a, H2O.port_a) annotation (Line(
            points={{-34,-16},{-34,16},{-44,16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(membrane2.port_a, HCO3.port_a) annotation (Line(
            points={{18,-18},{18,32},{24,32}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(membrane2.port_b, HCO3_E.port_a) annotation (Line(
            points={{18,-38},{18,-60},{26,-60}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(membrane3.port_b, Cl_E.port_a) annotation (Line(
            points={{72,-38},{72,-84},{78,-84}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(membrane3.port_a, Cl_P.port_a) annotation (Line(
            points={{72,-18},{72,12},{76,12}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(membrane.port_a, CO2_liquid.port_a) annotation (Line(
            points={{-68,-18},{-68,32},{-72,32}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(gasSolubility.gas_port, CO2_gas.port_a) annotation (Line(
            points={{-80,66},{-80,72}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HCO3.solution, blood_plasma.solution) annotation (Line(
            points={{40,22},{40,-16},{0,-16}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(HendersonHasselbalch1.products[2], H_E.port_a) annotation (Line(
            points={{4,-59.5},{14,-59.5},{14,-82},{48,-82}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>CO2 solution in water without any other acid-base buffers.</p>
<pre><b>plotExpression(apply(-log10(CarbonDioxideInWater.H3O.solute)),&nbsp;false,&nbsp;&QUOT;pH&QUOT;,&nbsp;1);</b></pre>
<p><br>Please note, that OH- (and CO3^-2) can be neglected from electroneutrality calculation, because of very small concentrations (in physiological pH) anyway. </p>
<p>And if SID&GT;0 then also H3O+ can be also neglected from electroneutrality, because only bicarbonate anions HCO3- (or CO3^-2) are needed there to balance the electroneutrality.</p>
<p><br>The partial pressure of CO2 in gas are input parameter. Outputs are an amount of free disolved CO2 in liquid and an amount of HCO3-.</p>
<p><br>The titration slope der(pH)/der(SID)=17.5 1/(mol/L) at pH=7.4 and pCO2=40 mmHg.</p>
<p><br>Molar heat of formation (aqueous):</p>
<p>CO2:        -413.5 kJ/mol  (gas: -393.5 kJ/mol )</p>
<p>H2O:        -285.8 kJ/mol</p>
<p>HCO3-:        -692.0 kJ/mol</p>
<p>CO3^-2:        -677.1 kJ/mol</p>
<p><br>Enthalphy of reaction H2O + CO2 &LT;-&GT; HCO3- + H+  :         7.3 kJ/mol</p>
<p>Enthalphy of reaction HCO3- &LT;-&GT; CO3^-2 + H+  :        14.9 kJ/mol</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=100, __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics),
          __Dymola_experimentSetupOutput);
      end CarbonDioxideInBlood;




    end AcidBase;

    package Hemoglobin "Hemoglobin blood gases binding"
      model Allosteric_Hemoglobin_MWC "Monod,Wyman,Changeux (1965)"
        extends Modelica.Icons.Example;

        constant Modelica.SIunits.AmountOfSubstance THb = 0.001
          "Total amount of hemoglobin";

        constant Modelica.SIunits.Temperature T=298.15 "Base Temperature";

        constant Real L(final unit="1")=7.0529*10^6
          "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
        constant Real c(final unit="1")=0.00431555
          "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
        constant Modelica.SIunits.Concentration KR=0.000671946
          "Oxygen dissociation coefficient on relaxed(R) hemoglobin subunit";

        constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 39.7
          "Amount of solution used for molarity to mole fraction conversion";

        constant Modelica.SIunits.Volume OneLiter = 0.001;
        constant Real KRx=KR*OneLiter/AmountOfSolutionIn1L
          "Mole fraction based KR";

      //Relative Gibbs formation energies of the substances in the system:
        constant Real RT = Modelica.Constants.R*T;
        constant Modelica.SIunits.MolarEnergy
          GO2aq=-RT*log(0.0013*0.018),
          GR0=0,                            GT0=GR0 -RT*log(L),
          GR1=GR0+GO2aq +RT*log(KRx/4),     GT1=GR1 -RT*log(c*L),
          GR2=GR1+GO2aq +RT*log(KRx/(3/2)), GT2=GR2 -RT*log(c^2*L),
          GR3=GR2+GO2aq +RT*log(KRx/(2/3)), GT3=GR3 -RT*log(c^3*L),
          GR4=GR3+GO2aq +RT*log(KRx*4),     GT4=GR4 -RT*log(c^4*L);

        parameter Modelica.SIunits.Time Tau = 1 "Slow down factor";

        Components.Substance oxygen_unbound(substanceData( DfG_25degC=GO2aq),
            amountOfSubstance_start(displayUnit="mol") = 1e-5)
          annotation (Placement(transformation(extent={{-56,-44},{-36,-24}})));

        Components.Substance T0(substanceData( DfG_25degC=GT0), amountOfSubstance_start=
              THb)
          annotation (Placement(transformation(extent={{34,78},{54,98}})));

        Components.Substance T1(substanceData( DfG_25degC=GT1),
            amountOfSubstance_start=THb*1e-4)
          annotation (Placement(transformation(extent={{34,36},{54,56}})));

        Components.Substance T2(substanceData( DfG_25degC=GT2),
            amountOfSubstance_start=THb*1e-8)
          annotation (Placement(transformation(extent={{34,-10},{54,10}})));

        Components.Substance R1(substanceData( DfG_25degC=GR1),
            amountOfSubstance_start=THb*1e-8)
          annotation (Placement(transformation(extent={{-20,36},{0,56}})));

        Components.Substance R2(substanceData( DfG_25degC=GR2),
            amountOfSubstance_start=THb*1e-10)
          annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

        Components.Substance T3(substanceData( DfG_25degC=GT3),
            amountOfSubstance_start=THb*1e-12)
          annotation (Placement(transformation(extent={{34,-54},{54,-34}})));

        Components.Substance R3(substanceData( DfG_25degC=GR3),
            amountOfSubstance_start=THb*1e-12)
          annotation (Placement(transformation(extent={{-20,-54},{0,-34}})));

        Components.Substance T4(substanceData( DfG_25degC=GT4),
            amountOfSubstance_start=THb*1e-17)
          annotation (Placement(transformation(extent={{34,-92},{54,-72}})));

        Components.Substance R4(substanceData( DfG_25degC=GR4),
            amountOfSubstance_start=THb*1e-14)
          annotation (Placement(transformation(extent={{-20,-92},{0,-72}})));

        Components.Substance R0(substanceData( DfG_25degC=GR0),
            amountOfSubstance_start=THb*1e-7)
          annotation (Placement(transformation(extent={{-20,78},{0,98}})));

        Components.Reaction quaternaryForm(Tau=Tau)
                                           annotation (Placement(transformation(extent={{4,78},{24,98}})));
        Components.Reaction oxyR1(nP=2, Tau=Tau)
                                        annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,64})));
        Components.Reaction oxyT1(nP=2, Tau=Tau)
                                        annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,64})));
        Components.Reaction oxyR2(nP=2, Tau=Tau)
                                        annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,22})));
        Components.Reaction oxyR3(nP=2, Tau=Tau)
                                        annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-24})));
        Components.Reaction oxyR4(nP=2, Tau=Tau)
                                        annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-66})));
        Components.Reaction oxyT2(nP=2, Tau=Tau)
                                        annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,22})));
        Components.Reaction oxyT3(nP=2, Tau=Tau)
                                        annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-24})));
        Components.Reaction oxyT4(nP=2, Tau=Tau)
                                        annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-66})));
        Components.Reaction quaternaryForm1(Tau=Tau)
                                            annotation (Placement(transformation(extent={{8,36},{28,56}})));
        Components.Reaction quaternaryForm2(Tau=Tau)
                                            annotation (Placement(transformation(extent={{8,-10},{28,10}})));
        Components.Reaction quaternaryForm3(Tau=Tau)
                                            annotation (Placement(transformation(extent={{8,-54},{28,-34}})));
        Components.Reaction quaternaryForm4(Tau=Tau)
                                            annotation (Placement(transformation(extent={{10,-92},{30,-72}})));

        Modelica.Blocks.Math.Sum oxygen_bound(k={1,1,2,2,3,3,4,4}, nin=8)
          annotation (Placement(transformation(extent={{72,-42},{82,-32}})));
        Modelica.Blocks.Math.Division sO2_ "hemoglobin oxygen saturation"
          annotation (Placement(transformation(extent={{86,-60},{96,-50}})));
        Modelica.Blocks.Math.Sum tHb(nin=10, k=4*ones(10))
          annotation (Placement(transformation(extent={{70,-80},{80,-70}})));

        Modelica.Blocks.Sources.Clock clock(offset=1000)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,62})));
        Sources.AirSubstance O2_in_air(
          TotalPressure(displayUnit="kPa") = 101325.0144354,
          substanceData = Chemical.Examples.Substances.Oxygen_gas,
          PartialPressure(displayUnit="kPa") = 1000,
          usePartialPressureInput=true) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,22})));

        Components.GasSolubility gasSolubility(useWaterCorrection=false, Tau=Tau/2)
          annotation (Placement(transformation(extent={{-94,-16},{-74,4}})));
        Components.Solution solution(amountOfSolution_start=
              AmountOfSolutionIn1L)
          annotation (Placement(transformation(extent={{-56,-102},{100,104}})));

      equation
       //  sO2 = (R1.amountOfSubstance + 2*R2.amountOfSubstance + 3*R3.amountOfSubstance + 4*R4.amountOfSubstance + T1.amountOfSubstance + 2*T2.amountOfSubstance + 3*T3.amountOfSubstance + 4*T4.amountOfSubstance)/(4*totalAmountOfHemoglobin);
      //   totalAmountOfRforms = R0.amountOfSubstance + R1.amountOfSubstance + R2.amountOfSubstance + R3.amountOfSubstance + R4.amountOfSubstance;
      //   totalAmountOfTforms = T0.amountOfSubstance + T1.amountOfSubstance + T2.amountOfSubstance + T3.amountOfSubstance + T4.amountOfSubstance;

      //   totalAmountOfHemoglobin*normalizedState[1] = totalAmountOfRforms + totalAmountOfTforms;

        connect(quaternaryForm.products[1],T0. port_a) annotation (Line(
            points={{24,88},{54,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR1.substrates[1],R1. port_a) annotation (Line(
            points={{-10,54},{-10,46},{0,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R1.port_a,oxyR2. products[1]) annotation (Line(
            points={{0,46},{0,32},{-10.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR2.substrates[1],R2. port_a) annotation (Line(
            points={{-10,12},{-10,0},{0,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.substrates[1],R3. port_a) annotation (Line(
            points={{-10,-34},{-10,-44},{0,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.products[1],R2. port_a) annotation (Line(
            points={{-10.5,-14},{-10.5,-7},{0,-7},{0,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R3.port_a,oxyR4. products[1]) annotation (Line(
            points={{0,-44},{0,-56},{-10.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR4.substrates[1],R4. port_a) annotation (Line(
            points={{-10,-76},{-10,-82},{0,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT1.products[1],T0. port_a) annotation (Line(
            points={{44.5,74},{44.5,88},{54,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT1.substrates[1],T1. port_a) annotation (Line(
            points={{44,54},{44,46},{54,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(T1.port_a,oxyT2. products[1]) annotation (Line(
            points={{54,46},{54,32},{44.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT3.substrates[1],T3. port_a) annotation (Line(
            points={{44,-34},{44,-44},{54,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(T3.port_a,oxyT4. products[1]) annotation (Line(
            points={{54,-44},{54,-56},{44.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT4.substrates[1],T4. port_a) annotation (Line(
            points={{44,-76},{44,-82},{54,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,quaternaryForm. substrates[1]) annotation (Line(
            points={{0,88},{4,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,oxyR1. products[1]) annotation (Line(
            points={{0,88},{0,74},{-10.5,74}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R1.port_a,quaternaryForm1. substrates[1]) annotation (Line(
            points={{0,46},{8,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm1.products[1],T1. port_a) annotation (Line(
            points={{28,46},{54,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R2.port_a,quaternaryForm2. substrates[1]) annotation (Line(
            points={{0,0},{8,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R3.port_a,quaternaryForm3. substrates[1]) annotation (Line(
            points={{0,-44},{8,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm3.products[1],T3. port_a) annotation (Line(
            points={{28,-44},{54,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R4.port_a,quaternaryForm4. substrates[1]) annotation (Line(
            points={{0,-82},{10,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm4.products[1],T4. port_a) annotation (Line(
            points={{30,-82},{54,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_bound.y,sO2_. u1) annotation (Line(
            points={{82.5,-37},{84,-37},{84,-52},{85,-52}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(sO2_.u2,tHb. y) annotation (Line(
            points={{85,-58},{84,-58},{84,-75},{80.5,-75}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(oxyR1.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,74},{-36,74},{-36,-34}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR2.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,32},{-36,32},{-36,-34}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,-14},{-36,-14},{-36,-34}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR4.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,-56},{-36,-56},{-36,-34}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT1.products[2])
                                            annotation (Line(
            points={{-36,-34},{-36,74},{43.5,74}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT2.products[2])
                                            annotation (Line(
            points={{-36,-34},{-36,32},{43.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT3.products[2])
                                            annotation (Line(
            points={{-36,-34},{-36,-14},{43.5,-14}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT4.products[2])
                                            annotation (Line(
            points={{-36,-34},{-36,-56},{43.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(O2_in_air.port_a, gasSolubility.gas_port) annotation (Line(
            points={{-84,12},{-84,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(gasSolubility.liquid_port, oxygen_unbound.port_a) annotation (Line(
            points={{-84,-16},{-84,-34},{-36,-34}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.solution, solution.solution) annotation (Line(
            points={{-52,-44},{-52,-102},{22,-102}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(R0.solution, solution.solution) annotation (Line(
            points={{-16,78},{-16,-102},{22,-102}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(T0.solution, solution.solution) annotation (Line(
            points={{38,78},{38,-102},{22,-102}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(R1.solution, solution.solution) annotation (Line(points={{-16,36},{-16,
                -102},{22,-102}}, smooth=Smooth.None));
        connect(T1.solution, solution.solution) annotation (Line(points={{38,36},{38,-102},
                {22,-102}}, smooth=Smooth.None));
        connect(R2.solution, solution.solution) annotation (Line(points={{-16,-10},{-16,
                -102},{22,-102}}, smooth=Smooth.None));
        connect(T3.solution, solution.solution) annotation (Line(points={{38,-54},{38,
                -102},{22,-102}}, smooth=Smooth.None));
        connect(R3.solution, solution.solution) annotation (Line(points={{-16,-54},{-16,
                -102},{22,-102}}, smooth=Smooth.None));
        connect(R4.solution, solution.solution) annotation (Line(points={{-16,-92},{-16,
                -102},{22,-102}}, smooth=Smooth.None));
        connect(T4.solution, solution.solution) annotation (Line(points={{38,-92},{38,
                -102},{22,-102}}, smooth=Smooth.None));
        connect(quaternaryForm2.products[1], T2.port_a) annotation (Line(
            points={{28,0},{54,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT2.substrates[1], T2.port_a) annotation (Line(
            points={{44,12},{44,0},{54,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(T2.port_a, oxyT3.products[1]) annotation (Line(
            points={{54,0},{54,-14},{44.5,-14}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(T2.solution, solution.solution) annotation (Line(points={{38,-10},{38,
                -102},{22,-102}}, smooth=Smooth.None));
        connect(R1.amountOfSubstance,oxygen_bound. u[1]) annotation (Line(
            points={{0,40},{64,40},{64,-37.875},{71,-37.875}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T1.amountOfSubstance,oxygen_bound. u[2]) annotation (Line(
            points={{54,40},{64,40},{64,-37.625},{71,-37.625}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(R2.amountOfSubstance,oxygen_bound. u[3]) annotation (Line(
            points={{0,-6},{64,-6},{64,-37.375},{71,-37.375}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(R3.amountOfSubstance,oxygen_bound. u[5]) annotation (Line(
            points={{0,-50},{64,-50},{64,-36.875},{71,-36.875}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T3.amountOfSubstance,oxygen_bound. u[6]) annotation (Line(
            points={{54,-50},{64,-50},{64,-36.625},{71,-36.625}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(R4.amountOfSubstance,oxygen_bound. u[7]) annotation (Line(
            points={{0,-88},{64,-88},{64,-36.375},{71,-36.375}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T4.amountOfSubstance,oxygen_bound. u[8]) annotation (Line(
            points={{54,-88},{64,-88},{64,-36.125},{71,-36.125}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T2.amountOfSubstance, oxygen_bound.u[4]) annotation (Line(
            points={{54,-6},{64,-6},{64,-37.125},{71,-37.125}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(R0.amountOfSubstance,tHb. u[1]) annotation (Line(
            points={{0,82},{64,82},{64,-75.9},{69,-75.9}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T0.amountOfSubstance,tHb. u[2]) annotation (Line(
            points={{54,82},{64,82},{64,-75.7},{69,-75.7}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(R1.amountOfSubstance,tHb. u[3]) annotation (Line(
            points={{0,40},{64,40},{64,-75.5},{69,-75.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T1.amountOfSubstance,tHb. u[4]) annotation (Line(
            points={{54,40},{64,40},{64,-75.3},{69,-75.3}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(R2.amountOfSubstance,tHb. u[5]) annotation (Line(
            points={{0,-6},{64,-6},{64,-75.1},{69,-75.1}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(R3.amountOfSubstance,tHb. u[7]) annotation (Line(
            points={{0,-50},{64,-50},{64,-74.7},{69,-74.7}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T3.amountOfSubstance,tHb. u[8]) annotation (Line(
            points={{54,-50},{64,-50},{64,-74.5},{69,-74.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(R4.amountOfSubstance,tHb. u[9]) annotation (Line(
            points={{0,-88},{64,-88},{64,-74.3},{69,-74.3}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T4.amountOfSubstance,tHb. u[10]) annotation (Line(
            points={{54,-88},{64,-88},{64,-74.1},{69,-74.1}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T2.amountOfSubstance, tHb.u[6]) annotation (Line(
            points={{54,-6},{64,-6},{64,-74.9},{69,-74.9}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(clock.y, O2_in_air.partialPressure) annotation (Line(
            points={{-84,51},{-84,32}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (          experiment(
            StopTime=15000,
            Tolerance=0.001,
            __Dymola_Algorithm="Dassl"),                  Documentation(info="<html>
<p>To understand the model is necessary to study the principles of MWC allosteric transitions first published by </p>
<p>[1] Monod,Wyman,Changeux (1965). &QUOT;On the nature of allosteric transitions: a plausible model.&QUOT; Journal of molecular biology 12(1): 88-118.</p>
<p><br>In short it is about binding oxygen to hemoglobin.</p>
<p>Oxgen are driven by its partial pressure using clock source - from very little pressure to pressure of 10kPa.</p>
<p>(Partial pressure of oxygen in air is the air pressure multiplied by the fraction of the oxygen in air.)</p>
<p>Hemoglobin was observed (by Perutz) in two structuraly different forms R and T.</p>
<p>These forms are represented by blocks T0..T4 and R0..R4, where the suffexed index means the number of oxygen bounded to the form.</p>
<p><br>In equilibrated model can be four chemical reactions removed and the results will be the same, but dynamics will change a lot. ;)</p>
<p>If you remove the quaternaryForm1,quaternaryForm2,quaternaryForm3,quaternaryForm4 then the model in equilibrium will be exactly the same as in MWC article.</p>
<p><br>Parameters was fitted to data of Severinghaus article from 1979. (For example at pO2=26mmHg is oxygen saturation sO2 = 48.27 &percnt;).</p>
</html>",   revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics),
          __Dymola_experimentSetupOutput);
      end Allosteric_Hemoglobin_MWC;

      model Allosteric_Hemoglobin2_MWC
        "Monod,Wyman,Changeux (1965) - The same allosteric hemoglobin model as Allosteric_Hemoglobin_MWC implemented by Speciation blocks"

       extends Modelica.Icons.Example;

        parameter Modelica.SIunits.MolarEnergy DfHT=10000
          "Enthalpy of formation of heme oxygenation in T hemoglobin form";
        parameter Modelica.SIunits.MolarEnergy DfHR=20000
          "Enthalpy of formation of heme oxygenation in R hemoglobin form";
        parameter Modelica.SIunits.MolarEnergy DfHL=-1000
          "Enthalpy of formation of reaction T->R as hemoglobin tetramer structure change";

        parameter Real L=7.0529*10^6
          "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
        parameter Modelica.SIunits.MoleFraction c=0.00431555
          "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
        parameter Modelica.SIunits.Concentration KR=0.000671946
          "oxygen dissociation on relaxed(R) hemoglobin subunit";
                                                                    //*7.875647668393782383419689119171e-5
                                                                  //10.500001495896 7.8756465463794e-05

        parameter Modelica.SIunits.Concentration KT=KR/c
          "oxygen dissociation on tensed(T) hemoglobin subunit";

        constant Real RT=Modelica.Constants.R*298.15;
        constant Modelica.SIunits.Volume OneLiter=0.001;

        parameter Modelica.SIunits.MoleFraction KRx = KR*OneLiter/55.508;
        parameter Modelica.SIunits.MoleFraction KTx = KT*OneLiter/55.508;

        parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508);
        parameter Modelica.SIunits.ChemicalPotential DfG_uR = 0;
        parameter Modelica.SIunits.ChemicalPotential DfG_uRO2 = DfG_uR + DfG_O2 + RT * log(KRx);
        parameter Modelica.SIunits.ChemicalPotential DfG_uT = 0;
        parameter Modelica.SIunits.ChemicalPotential DfG_uTO2 = DfG_uT + DfG_O2 + RT * log(KTx);
        parameter Modelica.SIunits.ChemicalPotential DfG_tT = 0;
        parameter Modelica.SIunits.ChemicalPotential DfG_tR = DfG_tT + RT * log(L);

        parameter Modelica.SIunits.AmountOfSubstance totalAmountOfHemoglobin=1;

        parameter Modelica.SIunits.Time Tau = 0.1 "Slow down factor";

        Components.Reaction quaternaryForm(Tau=Tau)
          annotation (Placement(transformation(extent={{0,-68},{20,-48}})));
        Components.Speciation R0_in_R(NumberOfSubunits=4, substanceData(DfG_25degC=DfG_tR),
          AmountOfSubstance_start=4e-11)
          annotation (Placement(transformation(extent={{-50,-68},{-30,-48}})));
        Components.Speciation T0_in_T(NumberOfSubunits=4, substanceData(DfG_25degC=DfG_tT),
          AmountOfSubstance_start=totalAmountOfHemoglobin)
          annotation (Placement(transformation(extent={{70,-66},{50,-46}})));
        Components.Substance OxyRHm[4](
          each amountOfSubstance_start=4e-19,
          substanceData( each DfH=-DfHL/4 - DfHR, each DfG_25degC=DfG_uRO2))
          "Oxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-96,-18},{-76,2}})));
        Components.Reaction oxygenation_R[4](each nP=2,each Tau=Tau)
          annotation (Placement(transformation(extent={{-68,-18},{-48,2}})));
        Components.Substance DeoxyRHm[4](
          each amountOfSubstance_start=4e-11,
          substanceData(each DfH=-DfHL/4, each DfG_25degC=DfG_uR))
          "Deoxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-40,-18},{-20,2}})));
        Components.Substance OxyTHm[4](
          substanceData(each DfH=-DfHT, each DfG_25degC=DfG_uTO2),
          each amountOfSubstance_start=1e-14)
          "Oxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{14,-18},{34,2}})));
        Components.Reaction oxygenation_T[4](each nP=2, each Tau=Tau)
          annotation (Placement(transformation(extent={{42,-18},{62,2}})));
        Components.Substance DeoxyTHm[4](
          substanceData(each DfH=0, each DfG_25degC=DfG_uT), each
            amountOfSubstance_start=totalAmountOfHemoglobin)
          "Deoxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{68,-18},{88,2}})));

        Components.Substance oxygen_unbound(        substanceData(DfG_25degC=DfG_O2),
            amountOfSubstance_start=1e-5)
          annotation (Placement(transformation(extent={{-2,6},{18,26}})));
        Modelica.Blocks.Sources.Clock clock(offset=1000)
          annotation (Placement(transformation(extent={{-40,74},{-20,94}})));
        Sources.AirSubstance        oxygen_in_air(
          usePartialPressureInput=true, substanceData=Chemical.Examples.Substances.Oxygen_gas)
                    annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={8,68})));
        Components.GasSolubility partialPressure1(Tau=Tau, useWaterCorrection=false)
                                                  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                origin={8,40})));
        Modelica.Blocks.Math.Sum sum1(nin=8, k=(1/4)*ones(8))
                                             annotation (Placement(transformation(
              extent={{-4,-4},{4,4}},
              rotation=270,
              origin={0,-36})));
        Components.Solution solution(                      ConstantTemperature=false)
          annotation (Placement(transformation(extent={{-100,-100},{100,56}})));

        Real sO2 "Hemoglobin oxygen saturation";
      equation
        sO2 = sum1.y/totalAmountOfHemoglobin;

        connect(OxyTHm.port_a, oxygenation_T.substrates[1])
                                                 annotation (Line(
            points={{34,-8},{42,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_T.products[1], DeoxyTHm.port_a)
                                               annotation (Line(
            points={{62,-8.5},{76,-8.5},{76,-8},{88,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(clock.y, oxygen_in_air.partialPressure) annotation (Line(
            points={{-19,84},{8,84},{8,78}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(OxyRHm.port_a, oxygenation_R.substrates[1]) annotation (Line(
            points={{-76,-8},{-68,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(DeoxyRHm.port_a, R0_in_R.subunits) annotation (Line(
            points={{-20,-8},{-40,-8},{-40,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_R.products[1], DeoxyRHm.port_a) annotation (Line(
            points={{-48,-8.5},{-34,-8.5},{-34,-8},{-20,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(T0_in_T.subunits, DeoxyTHm.port_a)   annotation (Line(
            points={{60,-46},{84,-46},{84,-8},{88,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(OxyRHm.amountOfSubstance, sum1.u[1:4]) annotation (Line(
            points={{-76,-14},{-70,-14},{-70,-28},{-0.1,-28},{-0.1,-31.2}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(OxyTHm.amountOfSubstance, sum1.u[5:8]) annotation (Line(
            points={{34,-14},{40,-14},{40,-28},{0.7,-28},{0.7,-31.2}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(oxygen_in_air.port_a, partialPressure1.gas_port) annotation (Line(
            points={{8,58},{8,50}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(partialPressure1.liquid_port, oxygen_unbound.port_a) annotation (Line(
            points={{8,30},{8,16},{18,16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(R0_in_R.port_a, quaternaryForm.substrates[1]) annotation (Line(
            points={{-30,-58},{0,-58}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm.products[1], T0_in_T.port_a) annotation (Line(
            points={{20,-58},{34,-58},{34,-56},{50,-56}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));

        for i in 1:4 loop
          connect(oxygenation_T[i].products[2], oxygen_unbound.port_a) annotation (Line(
            points={{62,-7.5},{78,-7.5},{78,16},{18,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
          connect(oxygenation_R[i].products[2], oxygen_unbound.port_a) annotation (Line(
            points={{-48,-7.5},{-12,-7.5},{-12,16},{18,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
          connect(OxyRHm[i].solution, solution.solution) annotation (Line(
            points={{-92,-18},{-92,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
          connect(DeoxyRHm[i].solution, solution.solution) annotation (Line(
            points={{-36,-18},{-92,-18},{-92,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
          connect(OxyTHm[i].solution, solution.solution) annotation (Line(
            points={{18,-18},{92,-18},{92,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
          connect(DeoxyTHm[i].solution, solution.solution) annotation (Line(
            points={{72,-18},{92,-18},{92,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        end for;

        connect(R0_in_R.solution, solution.solution) annotation (Line(
            points={{-46,-68},{-46,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(T0_in_T.solution, solution.solution) annotation (Line(
            points={{66,-66},{66,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(oxygen_unbound.solution, solution.solution) annotation (Line(points={{
                2,6},{92,6},{92,-100},{0,-100}}, smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics={Ellipse(
                extent={{2,12},{100,-36}},
                fillColor={255,181,181},
                fillPattern=FillPattern.Solid,
                pattern=LinePattern.None), Ellipse(
                extent={{-102,12},{-4,-36}},
                fillColor={255,181,181},
                fillPattern=FillPattern.Solid,
                pattern=LinePattern.None)}),
          experiment(StopTime=15000, __Dymola_Algorithm="Dassl"),
          Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Before silumation in &QUOT;Dymola 2014 FD01&QUOT; please chose &QUOT;Euler&QUOT; method!</p>
<p><br>To understand the model is necessary to study the principles of MWC allosteric transitions first published by </p>
<p>[1] Monod,Wyman,Changeux (1965). &QUOT;On the nature of allosteric transitions: a plausible model.&QUOT; Journal of molecular biology 12(1): 88-118.</p>
<p><br>In short it is about binding oxygen to hemoglobin.</p>
<p>Oxgen are driven by its partial pressure using clock source - from very little pressure to pressure of 10kPa.</p>
<p>(Partial pressure of oxygen in air is the air pressure multiplied by the fraction of the oxygen in air.)</p>
<p>Hemoglobin was observed (by Perutz) in two structuraly different forms R and T.</p>
<p>These forms are represented by blocks T0..T4 and R0..R4, where the suffexed index means the number of oxygen bounded to the form.</p>
<p><br>In equilibrated model can be four chemical reactions removed and the results will be the same, but dynamics will change a lot. ;)</p>
<p>If you remove the quaternaryForm1,quaternaryForm2,quaternaryForm3,quaternaryForm4 then the model in equilibrium will be exactly the same as in MWC article.</p>
<p><br>Parameters was fitted to data of Severinghaus article from 1979. (For example at pO2=26mmHg is oxygen saturation sO2 = 48.27 &percnt;).</p>
</html>"),__Dymola_experimentSetupOutput);
      end Allosteric_Hemoglobin2_MWC;

      model Hemoglobin_MKM_Specie "Part of model Hemoglobin_MKM_Adair"
         extends Interfaces.PartialSubstanceInSolution;

         parameter Modelica.SIunits.AmountOfSubstance amountOfSubstance_start = 6e-6
          "initial amount of substance";

      parameter Real[4] pKz
          "Dissociation coefficient of reaction z (Val1 amino terminal protonation)";
      parameter Real[4] pKc
          "Dissociation coefficient of reaction c (Val1 amino terminal carbamination)";
      parameter Real[4] pKh
          "Dissociation coefficient of reaction h (other Bohr protonation reactions of side chains)";

      protected
          constant Modelica.SIunits.Volume OneLiter = 0.001;
          constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 39.7
          "Amount of solution used for molarity to mole fraction conversion";

          parameter Real[4] Kz = {10,10,10,10}.^(-pKz)
          "Dissociation coefficient of reaction z (Val1 amino terminal protonation)";
          parameter Real[4] Kcx = ({10,10,10,10}.^(-pKc)) .* OneLiter./AmountOfSolutionIn1L
          "Dissociation coefficient of reaction c (Val1 amino terminal carbamination)";
          parameter Real[4] Kh = {10,10,10,10}.^(-pKh)
          "Dissociation coefficient of reaction h (other Bohr protonation reactions of side chains)";

          constant Real RT=Modelica.Constants.R*298.15;

          parameter Modelica.SIunits.MolarEnergy G_CO2 = Substances.CarbonDioxide_aqueous.DfG_25degC;
          parameter Modelica.SIunits.MolarEnergy G_A_H2[4] = {0,0,0,0};
          parameter Modelica.SIunits.MolarEnergy G_A_H3[4] = G_A_H2 + RT*log(Kz);
          parameter Modelica.SIunits.MolarEnergy G_A_HCOO[4] = G_A_H2 + G_CO2*ones(4) - RT*log(Kcx);
          parameter Modelica.SIunits.MolarEnergy G_AH_H2[4] = G_A_H2 - RT*log(Kh);
          parameter Modelica.SIunits.MolarEnergy G_AH_H3[4] = G_AH_H2 + RT*log(Kz);
          parameter Modelica.SIunits.MolarEnergy G_AH_HCOO[4] = G_AH_H2 + G_CO2*ones(4) - RT*log(Kcx);

          parameter Modelica.SIunits.MolarEnergy[4] dH_HbuANH2 = zeros(4)
          "Standard enthalpy of deprotonated and decarboxylated hemoglobin subunit";
      public
      parameter Modelica.SIunits.MolarEnergy[4] dHz
          "Enthalpy of reaction z (Val1 amino terminal protonation)";
      parameter Modelica.SIunits.MolarEnergy[4] dHc
          "Enthalpy of reaction c (Val1 amino terminal carbamination)";
      parameter Modelica.SIunits.MolarEnergy[4] dHh
          "Enthalpy of reaction h (other Bohr protonation reactions of side chains)";

          Components.Substance Hbu_A_NH3[4](
          substanceData(DfH=dH_HbuANH2 - dHz,  DfG_25degC=G_A_H3),
          each amountOfSubstance_start=1e-06)
          annotation (Placement(transformation(extent={{20,66},{40,86}})));
      Components.Substance Hbu_AH_NH3[4](
          each amountOfSubstance_start=1e-06,
          substanceData(DfH=dH_HbuANH2 - dHh - dHz,  DfG_25degC=G_AH_H3))
          annotation (Placement(transformation(extent={{-40,64},{-20,84}})));
      Components.Substance Hbu_A_NH2[4](
          substanceData(DfH=dH_HbuANH2,  DfG_25degC=G_A_H2), each
            amountOfSubstance_start=amountOfSubstance_start - 5e-6)
          annotation (Placement(transformation(extent={{20,-4},{40,16}})));
      Components.Substance Hbu_AH_NH2[4](
          each amountOfSubstance_start=1e-06,
          substanceData(DfH=dH_HbuANH2 - dHh,  DfG_25degC=G_AH_H2))
          annotation (Placement(transformation(extent={{-40,-8},{-20,12}})));
      Components.Substance Hbu_A_NHCOO[4](
          substanceData(DfH=dH_HbuANH2 + dHc,  DfG_25degC=G_A_HCOO),
          each amountOfSubstance_start=1e-06)
          annotation (Placement(transformation(extent={{20,-90},{40,-70}})));
      Components.Substance Hbu_AH_NHCOO[4](
          substanceData(DfH=dH_HbuANH2 + dHc,  DfG_25degC=G_AH_HCOO),
          each amountOfSubstance_start=1e-06)
          annotation (Placement(transformation(extent={{-40,-90},{-20,-70}})));
      Components.Reaction h2[4](
          each nS=1,
          each nP=2)
          annotation (Placement(transformation(extent={{-10,-4},{10,16}})));
          //K=fill(10, 4) .^ (-pKh .+ 3),
          //each TK=310.15,
          //dH=dHh
      Components.Reaction z1[4](
          each nP=2)
            annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,46})));
          //K=fill(10, 4) .^ (-pKz .+ 3),
          //dH=dHz,
          //each TK=310.15
      Components.Reaction z2[4](
          each nP=2)
            annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,40})));
          //K=fill(10, 4) .^ (-pKz .+ 3),
          //each TK=310.15,
          //dH=dHz
      Components.Reaction c1[4](
          each nS=2,
          each nP=2)
            annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,-36})));
          //K=fill(10, 4) .^ (-pKc .+ 3),
          //each TK=310.15,
          //dH=dHc
      Components.Reaction c2[4](
          each nS=2,
          each nP=2)
            annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-22,-40})));
          //K=fill(10, 4) .^ (-pKc .+ 3),
          //each TK=310.15,
          //dH=dHc
        Interfaces.SubstanceUsePort H "hydrogen ions"
          annotation (Placement(transformation(extent={{-110,70},{-90,90}})));
        Interfaces.SubstanceUsePort CO2
          annotation (Placement(transformation(extent={{-110,-70},{-90,-50}})));
        Components.Speciation Hb_tn(
          NumberOfSubunits=4, substanceData=substanceData,
          AmountOfSubstance_start=amountOfSubstance_start)
          annotation (Placement(transformation(extent={{66,-20},{86,0}})));
          Modelica.Blocks.Interfaces.RealOutput tHb_u(final unit="mol") annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={100,-60})));

              parameter Modelica.SIunits.Time Tau = 0.1 "Slow down factor";
      equation
       // x = Hb_tn.amountOfMacromolecule/AmountOfSolutionIn1L;
      connect(Hbu_AH_NH3.port_a, z2.substrates[1]) annotation (Line(
          points={{-20,74},{-20,50}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(Hbu_A_NH3.port_a, z1.substrates[1]) annotation (Line(
          points={{40,76},{40,56}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(z1.products[1], Hbu_A_NH2.port_a) annotation (Line(
          points={{39.5,36},{39.5,20},{40,20},{40,6}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(z2.products[1], Hbu_AH_NH2.port_a) annotation (Line(
          points={{-20.5,30},{-20.5,16},{-20,16},{-20,2}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(h2.substrates[1], Hbu_AH_NH2.port_a) annotation (Line(
          points={{-10,6},{-18,6},{-18,2},{-20,2}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(Hbu_A_NH2.port_a, c1.substrates[1]) annotation (Line(
          points={{40,6},{40,-26},{39.5,-26}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(Hbu_AH_NH2.port_a, c2.substrates[1]) annotation (Line(
          points={{-20,2},{-20,-30},{-22.5,-30}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(c1.products[1], Hbu_A_NHCOO.port_a) annotation (Line(
          points={{39.5,-46},{39.5,-80},{40,-80}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(c2.products[1], Hbu_AH_NHCOO.port_a) annotation (Line(
          points={{-22.5,-50},{-22.5,-66},{-20,-66},{-20,-80}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));

      connect(Hbu_A_NH2.port_a, h2.products[1]) annotation (Line(
          points={{40,6},{44,6},{44,5.5},{10,5.5}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));

        connect(Hb_tn.amountOfMacromolecule, tHb_u) annotation (Line(
            points={{78,-20},{78,-60},{100,-60}},
            color={0,0,127},
            smooth=Smooth.None));

        for i in 1:4 loop
          connect(z1[i].products[2], H) annotation (Line(
            points={{40.5,36},{40.5,26},{-56,26},{-56,80},{-100,80}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(z2[i].products[2], H) annotation (Line(
            points={{-19.5,30},{-19.5,26},{-56,26},{-56,80},{-100,80}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(c1[i].products[2], H) annotation (Line(
            points={{40.5,-46},{40.5,-54},{-56,-54},{-56,80},{-100,80}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(c2[i].products[2], H) annotation (Line(
            points={{-21.5,-50},{-21.5,-54},{-56,-54},{-56,80},{-100,80}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

          connect(H, h2[i].products[2]) annotation (Line(
            points={{-100,80},{-56,80},{-56,26},{14,26},{14,6.5},{10,6.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

          connect(CO2, c2[i].substrates[2]) annotation (Line(
            points={{-100,-60},{54,-60},{54,-16},{-21.5,-16},{-21.5,-30}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2, c1[i].substrates[2]) annotation (Line(
            points={{-100,-60},{54,-60},{54,-16},{40.5,-16},{40.5,-26}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(Hbu_A_NH2[i].solution, solution) annotation (Line(
            points={{24,-4},{24,-100},{-60,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(Hbu_A_NH3[i].solution, solution) annotation (Line(
            points={{24,66},{24,-100},{-60,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(Hbu_A_NHCOO[i].solution, solution) annotation (Line(
            points={{24,-90},{24,-100},{-60,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(Hbu_AH_NH3[i].solution, solution) annotation (Line(
            points={{-36,64},{-36,-100},{-60,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(Hbu_AH_NH2[i].solution, solution) annotation (Line(
            points={{-36,-8},{-36,-100},{-60,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(Hbu_AH_NHCOO[i].solution, solution) annotation (Line(
            points={{-36,-90},{-36,-100},{-60,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        end for;

        connect(Hb_tn.solution, solution) annotation (Line(
            points={{70,-20},{70,-100},{-60,-100}},
            color={0,0,0},
            smooth=Smooth.None));

        connect(Hbu_A_NH2.port_a, Hb_tn.subunits) annotation (Line(
            points={{40,6},{76,6},{76,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb_tn.port_a, port_a) annotation (Line(
            points={{86,-10},{100,-10},{100,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H, H) annotation (Line(
            points={{-100,80},{-100,80}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2014-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>

<p>[1] Mateják M, Kulhánek T, Matouaek S. Adair-Based Hemoglobin Equilibrium with Oxygen, Carbon Dioxide and Hydrogen Ion Activity. Scandinavian Journal of Clinical &AMP; Laboratory Investigation; 2015</p>

<p>[2] Bauer C, Schr&ouml;der E. Carbamino compounds of haemoglobin in human adult and foetal blood. The Journal of physiology 1972;227:457-71.</p>

<p>[3] Siggaard-Andersen O. Oxygen-Linked Hydrogen Ion Binding of Human Hemoglobin. Effects of Carbon Dioxide and 2, 3-Diphosphoglycerate I. Studies on Erythrolysate. Scandinavian Journal of Clinical &AMP; Laboratory Investigation 1971;27:351-60.</p>

</html>"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}),            graphics),
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics));
      end Hemoglobin_MKM_Specie;

      model Hemoglobin_MKM_Adair "Matejak,Kulhanek,Matousek (2014)"
        extends Modelica.Icons.Example;

        parameter Modelica.SIunits.Time Tau = 1e3 "Slow down factor";

        constant Real pKzD=7.73,pKcD=7.54,pKhD=4.52;
        constant Real pKzO=7.25,pKcO=8.35,pKhO=3.89;
        constant Modelica.SIunits.MolarEnergy dHzD=-51400;
        constant Modelica.SIunits.MolarEnergy dHzO=7700;
        constant Modelica.SIunits.MolarEnergy dHcD=59100;
        constant Modelica.SIunits.MolarEnergy dHcO=-41100;
        constant Modelica.SIunits.MolarEnergy dHhD=49000;
        constant Modelica.SIunits.MolarEnergy dHhO=-105000;
        constant Modelica.SIunits.MolarEnergy dHo=50000;
        constant Modelica.SIunits.MolarEnergy dH_HbuDANH2=0;
        // dHhD=0, dHhO=-104000, dHo=12700, dH_HbuDANH2=0;                           // dHhD=48600, dHhO=-104000, dHo=50000, dH_HbuDANH2=0;

        constant Modelica.SIunits.Volume OneLiter = 0.001;
        constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 39.7
          "Amount of solution used for molarity to mole fraction conversion";

        parameter Real K1x = 0.0121 *OneLiter/AmountOfSolutionIn1L;
        parameter Real K2x = 0.0117 *OneLiter/AmountOfSolutionIn1L;
        parameter Real K3x = 0.0871 *OneLiter/AmountOfSolutionIn1L;
        parameter Real K4x = 0.000386 *OneLiter/AmountOfSolutionIn1L;

      protected
         constant Real RT = Modelica.Constants.R*298.15;
        parameter Real GO2 = Chemical.Examples.Substances.Oxygen_aqueous.DfG_25degC;
        parameter Real G0 = 0;
        parameter Real G1 = G0 + GO2 + RT*log(K1x);
        parameter Real G2 = G1 + GO2 + RT*log(K2x);
        parameter Real G3 = G2 + GO2 + RT*log(K3x);
        parameter Real G4 = G3 + GO2 + RT*log(K4x);

      public
        Components.Reaction K1(
          nS=1,
          nP=2,
          Tau=Tau)
                annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-44,68})));
        Components.Reaction K2(
          nS=1,
          nP=2,
          Tau=Tau)
                annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-46,28})));
        Components.Reaction K3(
          nS=1,
          nP=2,
          Tau=Tau)
                annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-48,-18})));
        Components.Reaction K4(
          nS=1,
          nP=2,
          Tau=Tau)
                annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-50,-60})));
        Hemoglobin_MKM_Specie Hb0(
          substanceData(DfH = dH_HbuDANH2, DfG_25degC=G0),
          pKz=fill(pKzD, 4),
          pKc=fill(pKcD, 4),
          pKh=fill(pKhD, 4),
          dHz(displayUnit="kJ/mol") = fill(dHzD, 4),
          dHc(displayUnit="kJ/mol") = fill(dHcD, 4),
          dHh(displayUnit="kJ/mol") = fill(dHhD, 4),
          amountOfSubstance_start(displayUnit="mmol") = 0.00055,
          Tau=Tau)
          annotation (Placement(transformation(extent={{-4,78},{-24,98}})));
        Hemoglobin_MKM_Specie Hb1(
          substanceData(DfH = dH_HbuDANH2 + dHo, DfG_25degC=G1),
          pKz=cat(  1,
                    fill(pKzD, 3),
                    fill(pKzO, 1)),
          pKc=cat(  1,
                    fill(pKcD, 3),
                    fill(pKcO, 1)),
          pKh=cat(  1,
                    fill(pKhD, 3),
                    fill(pKhO, 1)),
          dHz(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dHzD, 3),
                  fill(dHzO, 1)),
          dHc(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dHcD, 3),
                  fill(dHcO, 1)),
          dHh(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dHhD, 3),
                  fill(dHhO, 1)),
          amountOfSubstance_start(displayUnit="mmol") = 0.00025,
          Tau=Tau)
          annotation (Placement(transformation(extent={{-4,40},{-24,60}})));

        Hemoglobin_MKM_Specie Hb2(
          substanceData(DfH = dH_HbuDANH2 + 2*dHo, DfG_25degC=G2),
          pKz=cat(  1,
                    fill(pKzD, 2),
                    fill(pKzO, 2)),
          pKc=cat(  1,
                    fill(pKcD, 2),
                    fill(pKcO, 2)),
          pKh=cat(  1,
                    fill(pKhD, 2),
                    fill(pKhO, 2)),
          dHz(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dHzD, 2),
                  fill(dHzO, 2)),
          dHc(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dHcD, 2),
                  fill(dHcO, 2)),
          dHh(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dHhD, 2),
                  fill(dHhO, 2)),
          amountOfSubstance_start(displayUnit="mmol") = 0.000106,
          Tau=Tau)
          annotation (Placement(transformation(extent={{-4,0},{-24,20}})));

        Hemoglobin_MKM_Specie Hb3(
          substanceData(DfH = dH_HbuDANH2 + 3*dHo, DfG_25degC=G3),
          pKz=cat(  1,
                    fill(pKzD, 1),
                    fill(pKzO, 3)),
          pKc=cat(  1,
                    fill(pKcD, 1),
                    fill(pKcO, 3)),
          pKh=cat(  1,
                    fill(pKhD, 1),
                    fill(pKhO, 3)),
          dHz(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dHzD, 1),
                  fill(dHzO, 3)),
          dHc(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dHcD, 1),
                  fill(dHcO, 3)),
          dHh(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dHhD, 1),
                  fill(dHhO, 3)),
          amountOfSubstance_start(displayUnit="mmol") = 8e-06,
          Tau=Tau)
          annotation (Placement(transformation(extent={{-4,-44},{-24,-24}})));
        Hemoglobin_MKM_Specie Hb4(
          substanceData(DfH = dH_HbuDANH2 + 4*dHo, DfG_25degC=G4),
          pKz=fill(pKzO, 4),
          pKc=fill(pKcO, 4),
          pKh=fill(pKhO, 4),
          dHz(displayUnit="kJ/mol") = fill(dHzO, 4),
          dHc(displayUnit="kJ/mol") = fill(dHcO, 4),
          dHh(displayUnit="kJ/mol") = fill(dHhO, 4),
          amountOfSubstance_start(displayUnit="mmol") = 9.6e-05,
          Tau=Tau)
          annotation (Placement(transformation(extent={{-4,-88},{-24,-68}})));
        Sources.AirSubstance        CO2(
                             substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
            PartialPressure(displayUnit="mmHg") = 5332.8954966)
          annotation (Placement(transformation(extent={{98,82},{78,102}})));
        Sources.AmbientMoleFraction      pH(substanceData=Chemical.Examples.Substances.Proton_aqueous,
            MoleFraction=10^(-7.4))                annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={40,84})));
        Modelica.Blocks.Math.Sum sO2(nin=4, k={4/4,3/4,2/4,1/4})
          annotation (Placement(transformation(extent={{32,-88},{52,-68}})));
        Components.Substance oxygen_unbound(
                     substanceData=Chemical.Examples.Substances.Oxygen_aqueous,
            amountOfSubstance_start=1e-05)
          annotation (Placement(transformation(extent={{-94,-28},{-74,-8}})));
        Modelica.Blocks.Sources.Clock clock(offset=1000)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,70})));
        Sources.AirSubstance        oxygen_in_air(
          PartialPressure(displayUnit="Pa") = 10,
          usePartialPressureInput=true,
          substanceData=Chemical.Examples.Substances.Oxygen_gas)
                    annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,34})));
        Components.GasSolubility partialPressure1(useWaterCorrection=true, Tau=Tau)
                                                  annotation (Placement(transformation(extent={{-10,-10},{10,
                  10}}, origin={-84,6})));
        Components.GasSolubility partialPressure2(Tau=Tau)
                                                  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                origin={68,62})));
        Components.Substance CO2_unbound(amountOfSubstance_start=
              0.0012, substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous)
          annotation (Placement(transformation(extent={{58,30},{78,50}})));
        Components.Solution solution(ConstantTemperature=false)
          annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
      equation
        connect(oxygen_unbound.port_a, K2.products[1]) annotation (Line(
            points={{-74,-18},{-62,-18},{-62,42},{-46,42},{-46,38},{-46.5,38}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, K3.products[1]) annotation (Line(
            points={{-74,-18},{-62,-18},{-62,0},{-48.5,0},{-48.5,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(K1.products[1], oxygen_unbound.port_a) annotation (Line(
            points={{-44.5,78},{-44.5,80},{-62,80},{-62,-18},{-74,-18}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, K4.products[1]) annotation (Line(
            points={{-74,-18},{-62,-18},{-62,-44},{-50.5,-44},{-50.5,-50}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(CO2_unbound.port_a, Hb0.CO2) annotation (Line(
            points={{78,40},{4,40},{4,82},{-4,82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb1.H, Hb0.H) annotation (Line(
            points={{-4,58},{10,58},{10,96},{-4,96}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb3.H, Hb0.H) annotation (Line(
            points={{-4,-26},{10,-26},{10,96},{-4,96}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb4.H, Hb0.H) annotation (Line(
            points={{-4,-70},{10,-70},{10,96},{-4,96}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb2.H, Hb0.H) annotation (Line(
            points={{-4,18},{10,18},{10,96},{-4,96}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb1.CO2, Hb0.CO2) annotation (Line(
            points={{-4,44},{4,44},{4,82},{-4,82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb2.CO2, Hb0.CO2) annotation (Line(
            points={{-4,4},{4,4},{4,82},{-4,82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb3.CO2, Hb0.CO2) annotation (Line(
            points={{-4,-40},{4,-40},{4,82},{-4,82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb4.CO2, Hb0.CO2) annotation (Line(
            points={{-4,-84},{4,-84},{4,82},{-4,82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb0.port_a, K1.products[2]) annotation (Line(
            points={{-24,88},{-43.5,88},{-43.5,78}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb1.port_a, K1.substrates[1]) annotation (Line(
            points={{-24,50},{-44,50},{-44,58}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb1.port_a, K2.products[2]) annotation (Line(
            points={{-24,50},{-44,50},{-44,38},{-45.5,38}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb2.port_a, K2.substrates[1]) annotation (Line(
            points={{-24,10},{-46,10},{-46,18}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb2.port_a, K3.products[2]) annotation (Line(
            points={{-24,10},{-46,10},{-46,-8},{-47.5,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb3.port_a, K3.substrates[1]) annotation (Line(
            points={{-24,-34},{-48,-34},{-48,-28}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb3.port_a, K4.products[2]) annotation (Line(
            points={{-24,-34},{-48,-34},{-48,-50},{-49.5,-50}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb4.port_a, K4.substrates[1]) annotation (Line(
            points={{-24,-78},{-50,-78},{-50,-70}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(Hb1.tHb_u, sO2.u[4]) annotation (Line(
            points={{-24,44},{-28,44},{-28,34},{16,34},{16,-76.5},{30,-76.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb2.tHb_u, sO2.u[3]) annotation (Line(
            points={{-24,4},{-32,4},{-32,-4},{18,-4},{18,-77.5},{30,-77.5}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(Hb3.tHb_u, sO2.u[2]) annotation (Line(
            points={{-24,-40},{-32,-40},{-32,-48},{20,-48},{20,-78.5},{30,-78.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb4.tHb_u, sO2.u[1]) annotation (Line(
            points={{-24,-84},{-32,-84},{-32,-96},{22,-96},{22,-79.5},{30,-79.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(clock.y, oxygen_in_air.partialPressure) annotation (Line(
            points={{-84,59},{-84,44}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb0.solution, solution.solution) annotation (Line(
            points={{-8,78},{-8,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(Hb1.solution, solution.solution) annotation (Line(points={{-8,40},{-8,
                -100},{0,-100}},      smooth=Smooth.None));
        connect(Hb2.solution, solution.solution) annotation (Line(points={{-8,0},{-8,-100},
                {0,-100}},            smooth=Smooth.None));
        connect(Hb3.solution, solution.solution) annotation (Line(points={{-8,-44},{-8,
                -100},{0,-100}},      smooth=Smooth.None));
        connect(Hb4.solution, solution.solution) annotation (Line(points={{-8,-88},{-8,
                -100},{0,-100}},      smooth=Smooth.None));
        connect(CO2_unbound.solution, solution.solution) annotation (Line(
            points={{62,30},{62,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(oxygen_unbound.solution, solution.solution) annotation (Line(
            points={{-90,-28},{-90,-100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(Hb0.H, pH.port_a) annotation (Line(
            points={{-4,96},{10,96},{10,84},{30,84}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(partialPressure2.gas_port, CO2.port_a) annotation (Line(
            points={{68,72},{68,92},{78,92}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_in_air.port_a, partialPressure1.gas_port) annotation (Line(
            points={{-84,24},{-84,16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(partialPressure1.liquid_port, oxygen_unbound.port_a) annotation (Line(
            points={{-84,-4},{-84,-18},{-74,-18}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(partialPressure2.liquid_port, CO2_unbound.port_a) annotation (Line(
            points={{68,52},{68,40},{78,40}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        annotation (          experiment(
            StopTime=15000,
            Tolerance=0.001,
            __Dymola_Algorithm="Dassl"), Documentation(info="<html>
<p>Before silumation in &QUOT;Dymola 2014 FD01&QUOT; please set environment variable &QUOT;<code><b>Advanced.Define.NonLinearIterations&nbsp;=&nbsp;3&QUOT;</b></code> and chose &QUOT;Euler&QUOT; method!</p>

<p>[1] Mateják M, Kulhánek T, Matouaek S. Adair-Based Hemoglobin Equilibrium with Oxygen, Carbon Dioxide and Hydrogen Ion Activity. Scandinavian Journal of Clinical &AMP; Laboratory Investigation; 2015</p>

<p>[2] Bauer C, Schr&ouml;der E. Carbamino compounds of haemoglobin in human adult and foetal blood. The Journal of physiology 1972;227:457-71.</p>

<p>[3] Siggaard-Andersen O. Oxygen-Linked Hydrogen Ion Binding of Human Hemoglobin. Effects of Carbon Dioxide and 2, 3-Diphosphoglycerate I. Studies on Erythrolysate. Scandinavian Journal of Clinical &AMP; Laboratory Investigation 1971;27:351-60.</p>

<p>[4] Severinghaus JW. Simple, accurate equations for human blood O2 dissociation computations. Journal of Applied Physiology 1979;46:599-602.</p>
</html>", revisions="<html>
<p><i>2014-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics),
          __Dymola_experimentSetupOutput);
      end Hemoglobin_MKM_Adair;



    end Hemoglobin;

  end Examples;

  package Components
    extends Modelica.Icons.Package;
    model Solution
      "Chemical solution as homogenous mixture of the substances at constant pressure"

      extends Interfaces.PartialSolution;

      parameter Modelica.SIunits.Pressure ConstantPressure=101325
        "Constant pressure of the solution";

      parameter Boolean ElectricGround = true
        "Is the solution electric potential equal to zero during simulation?";

      parameter Boolean ConstantTemperature = true
        "Has the solution constant temperature during simulation?";

    equation
      //isobaric condition
      solution.p = ConstantPressure;

      if  ElectricGround then
        //Solution connected to ground has zero voltage. However, electric current from the solution can varies.
        solution.v = 0;
      else
        //Electrically isolated solution has not any electric current from/to the solution. However, electric potential can varies.
        solution.i = 0;
      end if;

      if ConstantTemperature then
        //Ideal thermal exchange between environment and solution to reach constant temperature
        //0 = der(solution.T) = der((H-G)/S) = 0 = ((dH-dG)*S-(H-G)*dS)/(S*S)
        heatFromEnvironment = solution.dG - solution.dH + (freeEnthalpy-solution.G)*solution.dS/solution.S;
      else
        //Thermally isolated without any thermal exchange with environment
        heatFromEnvironment = 0;
      end if;

                                                                                                          annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false, initialScale=1, extent={{-100,-100},{100,100}}),
            graphics={
            Rectangle(
              extent={{-100,80},{100,-100}},
              lineColor={127,0,127},
              radius=10),
            Line(
              points={{-100,100},{-100,72},{-100,72}},
              color={127,0,127},
              smooth=Smooth.None),
            Line(
              points={{100,100},{100,80},{100,70}},
              color={127,0,127},
              smooth=Smooth.None),
                      Text(
              extent={{-88,-90},{92,-98}},
              lineColor={0,0,255},
              textString="%name",
              horizontalAlignment=TextAlignment.Left)}),
        Documentation(revisions="<html>
<p>2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<h4>amountOfSubstances = &int; MolarFlows</h4>
<h4>electricCharge = &int; ElectricCurrents</h4>
<h4>freeEnthalpy = &int; EnthalpyChanges</h4>
<h4>freeEntropy = &int; EntropyChanges</h4>
<h4>freeGibbsEnergy = &int; GibbsEnergyChanges</h4>
<h4>electricEnergy = &int; ElectricPowers</h4>
<p>Integration of all substances together into one homogenous mixture - the solution.</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end Solution;

    model Substance "Substance in solution"
      extends Icons.Substance;
      extends Interfaces.PartialSubstanceInSolution;

      parameter Modelica.SIunits.AmountOfSubstance amountOfSubstance_start=1e-8
        "Initial amount of the substance in compartment";

       Modelica.Blocks.Interfaces.RealOutput amountOfSubstance(start=amountOfSubstance_start, stateSelect=StateSelect.avoid, final unit="mol")
        "Current amount of the substance" annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,-60}),  iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,-60})));

      Real log10n(stateSelect=StateSelect.prefer)
        "Decadic logarithm of the amount of the substance in solution";
    protected
      constant Real InvLog_10=1/log(10);
    initial equation
      amountOfSubstance=amountOfSubstance_start;
    equation
      //der(amountOfSubstance)=port_a.q;
      //                                 log10n=log10(amountOfSubstance);
      //<- This is mathematically the same as two following lines. However, the differential solvers can handle the log10n much better. :-)
      der(log10n)=(InvLog_10)*(port_a.q/amountOfSubstance); amountOfSubstance = 10^log10n;

      //mole fraction (an analogy of molar concentration or molality)
      //if you select the amount of solution per one kilogram of solvent then the values of amountOfSubstance will be the same as molality
      //if you select the amount of solution in one liter of solution then the values of amountOfSubstance will be the same as molarity
      x = amountOfSubstance/solution.n;

      //local changes of the solution
      solution.dH = molarEnthalpy*port_a.q;
      solution.dS = molarEntropy*port_a.q;
      solution.dG = port_a.u*port_a.q;
      solution.dn = port_a.q;
      solution.i = Modelica.Constants.F * z * port_a.q;
      solution.dI = (1/2) * port_a.q * z^2;
      solution.dV = molarVolume * port_a.q;

                                                                                                         annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={Text(
              extent={{-80,90},{280,130}},
              lineColor={0,0,255},
              textString="%name")}),
        Documentation(revisions="<html>
<p>2009-2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<h4>n = x &middot; n(solution) = &int; MolarFlow</h4>
<p>where n is amount of the substance and x is mole fraction.</p>
<p>The main class from &ldquo;Chemical&rdquo; package is called &QUOT;Substance&QUOT;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &QUOT;n&QUOT; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &QUOT;n(solution)&QUOT; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solution. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
<p><br>The recalculation between mole fraction, molarity and molality can be written as follows:</p>
<p>x = n/n(solution) = b * m(solvent)/n(solution) = c * V(solution)/n(solution)</p>
<p>where m(solvent) is mass of solvent, V(solution) is volume of solution, b=n/m(solvent) is molality of the substance, c=n/V(solution) is molarity of the substance.</p>
<p>If the amount of solution is selected to the number of total solution moles per one kilogram of solvent then the values of x will be the same as molality.</p>
<p>If the amount of solution is selected to the number of total solution moles in one liter of solution then the values of x will be the same as molarity.</p>
<p><br><br>Definition of electro-chemical potential:</p>
<h4>u(x,T,V) = u&deg;(T) + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. relative eletric potential of the solution</p>
<p>z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
<p>R .. gas constant</p>
<p>F .. Faraday constant</p>
<p>gamma .. activity coefficient</p>
<p>u&deg;(T) .. chemical potential of pure substance</p>
<p>DfG(T) .. free Gibbs energy of formation of the substance at current temperature T. </p>
<p>DfH .. free enthalpy of formation of the substance</p>
<p>DfS .. free entropy of formation of the substance </p>
<p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end Substance;

    model Reaction "Chemical Reaction"

      parameter Integer nS=1 "Number of substrates types"
        annotation ( HideResult=true);

      parameter Modelica.SIunits.StoichiometricNumber s[nS]=ones(nS)
        "Stoichiometric reaction coefficient for substrates"
        annotation (HideResult=true);

      parameter Integer nP=1 "Number of products types"
        annotation ( HideResult=true);

      parameter Modelica.SIunits.StoichiometricNumber p[nP]=ones(nP)
        "Stoichiometric reaction coefficients for products"
        annotation (HideResult=true);

      Modelica.SIunits.MolarFlowRate rr(start=0) "Reaction molar flow rate";

      Interfaces.SubstanceUsePort products[nP] "Products"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      Interfaces.SubstanceUsePort substrates[nS] "Substrates"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

      //solution properties:
      parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution = 1
        "Amount of all particles in the reacting solution";

    /*  //for debugging:
  Real DissociationConstant "Dissociation constant as ratio of mole fractions";

  Modelica.SIunits.MolarEnergy DrH 
    "Standard Enthalpy Change of reaction (negative=exothermic)";

  Modelica.SIunits.Power lossHeat "Comsumed heat by the reaction";

//  Modelica.SIunits.ElectricPotential StandardNernstPotential
//    "Standard electric potential of half-cell rection";
*/
      parameter Modelica.SIunits.MolarEnergy ActivationEnergy(displayUnit="kJ/mol")=10000
        "To determine the reaction rate";

      parameter Modelica.SIunits.Time Tau = 1
        "Time constant for other scaling of the reaction rate";
    equation
      //the main equation
      rr = (AmountOfSolution/(Tau*ActivationEnergy))*((p * products.u) - (s * substrates.u));

      //reaction molar rates
      rr*s = -substrates.q;
      rr*p = products.q;

      //properties of solution
    //  amountOfSolution = substrates[1].amountOfSolution;

    /*  //for debugging olny:

    //evaluation of dissociation constant from the Gibbs energy of the reaction
    DissociationConstant = (product(substrates.activityCoefficient .^ s) / product(products.activityCoefficient .^ p)) *
    exp( -(1/(Modelica.Constants.R*substrates[1].temperature)) * (p * products.uPure - s * substrates.uPure));

  //  0 = p * products.u0 - s * substrates.u0 + (p*products.z - s*substrates.z) * Modelica.Constants.F * StandardNernstPotential;

    //molar enthalpy of reaction
    DrH = sum(p.*products.molarEnthalpy) - sum(s.*substrates.molarEnthalpy);

    //consumed heat by reaction
    lossHeat = DrH*rr; //DrH<0 => Exothermic => lossHeat>0, Endothermic otherwise
    */

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                            graphics={
            Rectangle(
              extent={{-100,-30},{100,30}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-128,-66},{142,-40}},
              textString="%name",
              lineColor={0,0,255}),
            Polygon(
              points={{-60,6},{-60,4},{54,4},{54,4},{18,14},{18,6},{-60,6}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{54,-8},{54,-6},{-60,-6},{-60,-6},{-24,-16},{-24,-8},{54,-8}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid)}),
        Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &LT;-&GT; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</b></sub> </p>
<p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
<p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
<td><p>dissociation constant</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>G= &Delta;<sub>r</sub>H - T&middot;&Delta;<sub>r</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K) </p></td>
<td><p>molar Gibb&apos;s energy of the reaction</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>H = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>H<sub>i</sub>) </p></td>
<td><p>molar enthalpy of the reaction</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>S = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>S<sub>i</sub>) = <a href=\"modelica://Modelica.Constants\">k</a>&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(&Delta;<sub>r</sub>&omega;) </p></td>
<td><p>molar entropy of the reaction</p></td>
</tr>
</table>
<p><br><h4><span style=\"color:#008000\">Notations</span></h4></p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>A<sub>i</sub></p></td>
<td><p>i-th substance</p></td>
</tr>
<tr>
<td><p>v<sub>i</sub></p></td>
<td><p>stochiometric coefficients of i-th substance</p></td>
</tr>
<tr>
<td><p>K</p></td>
<td><p>dissociation constant (activity based)</p></td>
</tr>
<tr>
<td><p>a(A<sub>i</sub>)=f<sub>i</sub>*x<sub>i</sub></p></td>
<td><p>activity of the substance A</p></td>
</tr>
<tr>
<td><p>f<sub>i</sub></p></td>
<td><p>activity coefficient of the substance A</p></td>
</tr>
<tr>
<td><p>x<sub>i</sub></p></td>
<td><p>mole fraction of the substance A</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>H<sub>i</sub></p></td>
<td><p>molar enthalpy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>G<sub>i</sub></p></td>
<td><p>molar Gibbs energy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>S<sub>i</sub></p></td>
<td><p>molar entropy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>&omega;</p></td>
<td><p>change of number of microstates of particles by reaction</p></td>
</tr>
<tr>
<td></td>
<td></td>
</tr>
</table>
</html>"));
    end Reaction;

    model Diffusion "Solute diffusion"
      extends Icons.Diffusion;
      extends Interfaces.OnePortParallel;

      parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution = 1
        "Amount of all particles in the diffusing solution";

      parameter Modelica.SIunits.MolarEnergy ActivationEnergy(displayUnit="kJ/mol")=10000
        "To determine the diffusion rate";

      parameter Modelica.SIunits.Time Tau = 1
        "Time constant for other scaling of the diffusion rate";

    equation
      port_b.q = (AmountOfSolution/(Tau*ActivationEnergy)) * (port_b.u - port_a.u);

       annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
    end Diffusion;

    model GasSolubility "Henry's law of gas solubility in liquid."

      extends Icons.GasSolubility;

      parameter Boolean useWaterCorrection = true
        "Are free Gibbs energy of aqueous formation shifted by 10 kJ/mol?"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));

      Interfaces.SubstanceUsePort gas_port "Gaseous solution"
        annotation (Placement(transformation(extent={{-10,90},{10,110}})));

      Interfaces.SubstanceUsePort liquid_port "Dissolved in liquid solution"
        annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
            iconTransformation(extent={{-10,-110},{10,-90}})));

            /*
  //for debugging
  Real kH(final unit="(mol/mol)/(Pa/Pa)", displayUnit="(mol/kg H2O)/bar at 25degC") 
    "Henry's law coefficient such as liquid-gas concentration ratio at 25degC";
  Modelica.SIunits.Temperature C(displayUnit="K") 
    "Henry's law temperature dependence coefficient";

  Modelica.SIunits.Power lossHeat "Comsumed heat by the reaction";

//  Modelica.SIunits.ElectricPotential StandardNernstPotential
//    "Standard electric potential";
*/
      parameter Modelica.SIunits.AmountOfSubstance AmountOfSubstanceInSurface = 1
        "The amount of liquid-gaseous molecules in the liquid-gas junction to the dissolution rate";

      parameter Modelica.SIunits.MolarEnergy ActivationEnergy(displayUnit="kJ/mol")=10000
        "To determine the diffusion rate";

      parameter Modelica.SIunits.Time Tau = 1
        "Time constant for other scaling of the diffusion rate";

    equation
      gas_port.q + liquid_port.q = 0;

      // the main equation
      liquid_port.q = (AmountOfSubstanceInSurface/(Tau*ActivationEnergy))*(liquid_port.u - gas_port.u - (if useWaterCorrection then Modelica.Constants.R*(298.15)*log(0.018) else 0));

      //for debugging olny:
    /*
    // evaluation of kH and C from enthalpies and Gibbs energies
    C=-(liquid_port.molarEnthalpy - gas_port.molarEnthalpy)/Modelica.Constants.R;

    -Modelica.Constants.R*liquid_port.temperature*
     log(kH*(liquid_port.activityCoefficient/gas_port.activityCoefficient)) =
     (liquid_port.uPure - gas_port.uPure)
     - (if useWaterCorrection then Modelica.Constants.R*liquid_port.temperature*log(0.018) else 0);

 //   0 = liquid_port.u0 - gas_port.u0 + (liquid_port.z - gas_port.z) * Modelica.Constants.F * StandardNernstPotential
 //    - (if useWaterCorrection then Modelica.Constants.R*liquid_port.temperature*log(0.018) else 0);

     lossHeat = (liquid_port.molarEnthalpy - gas_port.molarEnthalpy)*liquid_port.q; //negative = heat are comsumed when change from liquid to gas
     */

       annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Gaseuous substance dissolition in liquid (Henry&apos;s law, Raoult&apos;s law, Nernst dissolution in one). </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K<sub>H</sub> =x<sub>L</sub> / x<sub>g</sub>&nbsp;</p></td>
<td><p>Henry&apos;s coefficient</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>G = &Delta;<sub>sol</sub>H - T&middot;&Delta;<sub>sol</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K<sub>H</sub>&middot; (f<sub>L</sub> / f<sub>g</sub>)) </p></td>
<td><p>molar Gibb&apos;s energy of the dissolition</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>H = &Delta;<sub>f</sub>H<sub>L </sub>- &Delta;<sub>f</sub>H<sub>g</sub></p></td>
<td><p>molar enthalpy of the dissolition</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>S = &Delta;<sub>f</sub>S<sub>L</sub> - &Delta;<sub>f</sub>S<sub>g</sub> = <a href=\"modelica://Modelica.Constants\">k</a>&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(&Delta;<sub>sol</sub>&omega;) </p></td>
<td><p>molar entropy of the dissolition</p></td>
</tr>
</table>
<p><br><h4><span style=\"color:#008000\">Notations</span></h4></p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>x<sub>L</sub></p></td>
<td><p>mole fraction of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>x<sub>g</sub></p></td>
<td><p>mole fraction of the substance in the gas</p></td>
</tr>
<tr>
<td><p>f<sub>L</sub></p></td>
<td><p>activity coefficient of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>f<sub>g</sub></p></td>
<td><p>activity coefficient of the substance in the gas</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>H<sub>L</sub></p></td>
<td><p>molar enthalpy of formation of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>H<sub>g</sub></p></td>
<td><p>molar enthalpy of formation of the substance in the gas</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>S<sub>L</sub></p></td>
<td><p>molar entropy of formation of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>S<sub>g</sub></p></td>
<td><p>molar entropy of formation of the substance in the gas</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>G</p></td>
<td><p>molar Gibbs energy of dissolvation of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>&omega;</p></td>
<td><p>change of number of microstates of particles by dissolution</p></td>
</tr>
<tr>
<td></td>
<td></td>
</tr>
</table>
</html>"),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end GasSolubility;

    model Membrane
      "Passive transport of the substance through semipermeable membrane"
      extends Icons.Membrane;
      extends Interfaces.OnePortParallel;

      parameter Modelica.SIunits.AmountOfSubstance AmountOfMembraneChannels = 1
        "The amount of membrane channels to determine the rate of transport";

      parameter Modelica.SIunits.MolarEnergy ActivationEnergy(displayUnit="kJ/mol")=10000
        "To determine the transport rate";

      parameter Modelica.SIunits.Time Tau = 1
        "Time constant for other scaling of the transport rate";
     /*
  //for debugging
   parameter Modelica.SIunits.MolarVolume Vm=18.1367e-6 
    "Molar volume of the particle, defaultly set to water molar volume at 37degC";

   
   Modelica.SIunits.MoleFraction donnanRatio 
    "Donnan's ratios  as  x(inside)/x(outside)";

   Modelica.SIunits.OsmoticPressure opi 
    "Osmotic pressure of the substance on inner side of membrane";
   Modelica.SIunits.OsmoticPressure opo 
    "Osmotic pressure of the substance on outer side of membrane";

  Modelica.SIunits.ElectricPotential membranePotential 
    "Current potential on membrane";

  Modelica.SIunits.Power lossHeat "Comsumed heat by the reaction";

  Modelica.SIunits.ElectricPotential NernstPotential 
    "Nernst electric potential of the substance";
*/
    equation
      //the main equation
      port_a.q = (AmountOfMembraneChannels / (Tau*ActivationEnergy)) * (port_a.u - port_b.u);

      //for debuging only:
    /*
   //osmotic pressures of the substances
   opi = port_a.u ./ Vm;
   opo = port_b.u ./ Vm;

   //heat comsumption
   lossHeat = (port_a.molarEnthalpy - port_b.molarEnthalpy) * port_b.q;

   //evaluating Donnan's ratio on membrane
   Modelica.Constants.R * port_a.temperature * log(donnanRatio*(port_b.activityCoefficient/port_a.activityCoefficient))=
     - Modelica.Constants.F * membranePotential;

   //current potential on membrane
   membranePotential = port_a.electricPotential - port_b.electricPotential;

   //Nernst electric potential of the substance
   0 = (port_a.u - port_b.u) - (port_a.uPure - port_b.uPure) + port_a.z * Modelica.Constants.F * NernstPotential;
*/
      annotation ( Documentation(info="<html>
<p><u><b><font style=\"color: #008000; \">Filtration throught semipermeable membrane.</font></b></u></p>
<p>The penetrating particles are driven by electric and chemical gradient to reach Donnan&apos;s equilibrium.</p>
<p>If zero-flow Donnan&apos;s equilibrium is reached. </p>
</html>", revisions="<html>
<p><i>2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),     Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end Membrane;

    model Speciation
      "Quaternary form of macromolecule with independent subunits"
      extends Icons.Speciation;
      extends Interfaces.PartialSubstance;

      parameter Integer NumberOfSubunits=1
        "Number of independent subunits occuring in macromolecule";

      Interfaces.SubstanceUsePort subunits[NumberOfSubunits]
        "Subunits of macromolecule"
        annotation (Placement(transformation(extent={{-10,90},{10,110}})));

      Modelica.Blocks.Interfaces.RealOutput amountOfMacromolecule(final unit="mol")
        "Total amount of macromolecules including all selected forms of subunits"
                                                                                  annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-80}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={20,-100})));

      parameter Modelica.SIunits.AmountOfSubstance AmountOfSubstance_start=1e-8
        "Initial value of the total amount of the macromolecules. It must be the same as the total amount of each its subunit!";

       Real log10n(stateSelect=StateSelect.prefer)
        "Decadic logarithm of the amount of the substance in solution";

      Interfaces.SolutionPort solution annotation (Placement(transformation(extent={{-70,
                -110},{-50,-90}}),
            iconTransformation(extent={{-70,-110},{-50,-90}})));

    //  Real fractions[NumberOfSubunits]
    //    "Fractions of selected specific form of each subunit in macromolecule";
        parameter Modelica.SIunits.MolarEnergy ActivationEnergy(displayUnit="kJ/mol")=10000
        "To determine the speed of the equilibration";

        parameter Modelica.SIunits.Time Tau = 1
        "Time constant for other scaling of the speed of the reaction";
    protected
        Modelica.SIunits.MoleFraction xm
        "Mole fraction of all form of the macromolecule (in the conformation)";
        Modelica.SIunits.ChemicalPotential uEq
        "Chemical potential of the specific form of the macromolecule (in the conformation) at equilibrium";
    initial equation
      amountOfMacromolecule = AmountOfSubstance_start;
    equation
     // der(amountOfMacromolecule) = port_a.q;
     //                                        log10n=log10(amountOfMacromolecule);
      //<- This is mathematically the same as two following lines. However, the differential solvers can handle the log10n much better. :-)
      der(log10n)=port_a.q/(log(10)*amountOfMacromolecule);
      amountOfMacromolecule = 10^log10n;

      //change of macromolecule = change of its subunits
      subunits.q = -port_a.q * ones(NumberOfSubunits);

      port_a.q = (solution.n/(Tau*ActivationEnergy)) * (uEq - port_a.u);

      xm = amountOfMacromolecule/solution.n;
      uEq = substanceModel.u0(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength)
      + Modelica.Constants.R*temperature*log(xm) + sum(subunits.u - Modelica.Constants.R*temperature*log(xm)*ones(NumberOfSubunits));

    /*  //the amount of total macromolecule is the same as amount of each its selected subunit
  amountOfMacromolecule = amountOfSubunits[1];

  //chemical speciation
  x = (amountOfMacromolecule/solution.n)*product(fractions);

  fractions = if (amountOfMacromolecule < Modelica.Constants.eps) then zeros(NumberOfSubunits)
              else subunits.x ./ (amountOfSubunits/solution.n);
              */

      //global properties of the solution
      temperature = solution.T;
      pressure = solution.p;
      electricPotential = solution.v;
      amountOfSolution = solution.n;
      moleFractionBasedIonicStrength = solution.I;

      //changes of the solution, where all subunits are also connected
      solution.dH = molarEnthalpy * port_a.q; // 0 = subunits.molarEnthalpy * subunits.q;
      solution.dS = molarEntropy * port_a.q; // 0 = subunits.molarEntropy * subunits.q;
      solution.dG = port_a.u * port_a.q + subunits.u * subunits.q;
      solution.dn = port_a.q + sum(subunits.q);
      solution.i = 0; // 0 = Modelica.Constants.F * (port_a.z * port_a.q + subunits.z * subunits.q);
      solution.dI = 0; // 0 = (1/2) * port_a.q * port_a.z^2 + (1/2) * subunits.q * (subunits.z .^ 2);
      solution.dV = molarVolume * port_a.q; //0 = subunits.molarVolume * subunits.q;
      annotation (defaultComponentName="macromolecule",
        Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>Macromolecule speciation in chemical equilibrium</b> </p>
<p>The equilibrium of the conformation reactions of macromolecules can be simplified to the reactions of their selected electro-neutral forms of the selected conformation, because of the law of detailed balance.</p>
<p>The assumptions of this calculation are:</p>
<ol>
<li><i>Initial total concentrations of each subunit must be set to the total macromolecule concentration (for the selected conformation).</i></li>
<li><i>The charge, enthalpy of formation, entropy of formation and molar volume of each selected independent subunit form is zero. </i></li>
<li><i>Subunits are connected to the same solution as the macromolecule. </i></li>
</ol>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>x<sub>m</sub>&nbsp;</p></td>
<td><p>the probability of macromolecule(of the selected conformation)</p></td>
</tr>
<tr>
<td><p>f<sub>i</sub> = (x<sub>i</sub>/x<sub>m</sub>)</p></td>
<td><p>the probalitivy of selected independent subunits forms (of the macromolecule in the selected conformation)</p></td>
</tr>
<tr>
<td><p>x<sub>s </sub>= x<sub>m</sub>&middot; &Pi; f<sub>i</sub> = x<sub>m</sub>&middot; &Pi; (x<sub>i</sub>/x<sub>m</sub>)</p></td>
<td><p>the probability of the selected form of macromolecule (composed from selected subunits in the selected conformation)</p></td>
</tr>
<tr>
<td><p>u<sub>s </sub>= u<sub>s</sub>&deg; + R&middot;T&middot;ln(x<sub>m</sub>) + &sum; (u<sub>i</sub> - R&middot;T&middot;ln(x<sub>m</sub>))</p></td>
<td><p>final equation of the equilibrium of electro-chemical potential</p></td>
</tr>
</table>
<p><br><br><br><b><font style=\"color: #008000; \">Notations</font></b></p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>n<sub>T</sub></p></td>
<td><p>total amount of substances in the solution</p></td>
</tr>
<tr>
<td><p>n<sub>m</sub></p></td>
<td><p>total amount of the macromolecule (of the selected conformation) in the solution</p></td>
</tr>
<tr>
<td><p>n<sub>s</sub></p></td>
<td><p>amount of the specific form of the macromolecule (of the selected conformation) in the solution</p></td>
</tr>
<tr>
<td><p>n<sub>i</sub></p></td>
<td><p>amount of the specific form of the i-th macromolecule(of the selected conformation)&apos;s subunit in the solution</p></td>
</tr>
<tr>
<td><p>x<sub>m </sub>= n<sub>m </sub>/ n<sub>T</sub></p></td>
<td><p>mole fraction of macromolecule (of the selected conformation)</p></td>
</tr>
<tr>
<td><p>x<sub>s </sub>= n<sub>s </sub>/ n<sub>T</sub></p></td>
<td><p>mole fraction of the selected form of the whole macromolecule (composed from selected subunits in the selected conformation)</p></td>
</tr>
<tr>
<td><p>x<sub>i </sub>= n<sub>i </sub>/ n<sub>T</sub></p></td>
<td><p>mole fraction of i-th macromolecule(of the selected conformation)&apos;s subunit form</p></td>
</tr>
<tr>
<td><p>u<sub>s</sub>&deg;</p></td>
<td><p>base chemical potential of the selected form of the macromolecule (composed from selected subunits in the selected conformation)</p></td>
</tr>
<tr>
<td><p>u<sub>s </sub>= u<sub>s</sub>&deg; + R&middot;T&middot;ln(x<sub>s</sub>)</p></td>
<td><p>chemical potential of the selected form of the macromolecule (composed from selected subunits in the selected conformation)</p></td>
</tr>
<tr>
<td><p>u<sub>i</sub>&deg; = 0</p></td>
<td><p>base chemical potential of the specific form of the i-th macromolecule(of the selected conformation)&apos;s subunit in the solution</p></td>
</tr>
<tr>
<td><p>u<sub>i </sub>= R&middot;T&middot;ln(x<sub>i</sub>)</p></td>
<td><p>chemical potential of the specific form of the i-th macromolecule(of the selected conformation)&apos;s subunit in the solution</p></td>
</tr>
</table>
<p><br><br><br><br>For example: If the macromolecule M has four identical independent subunits and each subunit can occur in two form F1 and F2, then the probability of macromolecule form S composed only from four subunits in form F1 is P(S)=P(M)*P(F1)^4.</p>
</html>"),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics={                                                        Text(
              extent={{-22,-106},{220,-140}},
              lineColor={0,0,255},
              textString="%name")}),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end Speciation;

    model Stream "Flow of whole solution"
      extends Interfaces.OnePortParallel;
      extends Interfaces.ConditionalSolutionFlow;
      extends Interfaces.PartialSubstanceNoStorage;

    equation
      port_a.q = if (q>0) then q*x else q*x;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics={
            Rectangle(
              extent={{-100,-50},{100,50}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              rotation=360),
            Polygon(
              points={{-80,25},{80,0},{-80,-25},{-80,25}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              rotation=360),
            Text(
              extent={{-150,-20},{150,20}},
              textString="%name",
              lineColor={0,0,255},
              origin={2,-74},
              rotation=180)}),
        Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><h4><font color=\"#008000\">Bidirectional mass flow by concentration</font></h4></p>
<p>Possible field values: </p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0.1\"><tr>
<td></td>
<td><p align=\"center\"><h4>forward flow</h4></p></td>
<td><p align=\"center\"><h4>backward flow</h4></p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>solutionFlow</h4></p></td>
<td><p align=\"center\">&GT;=0</p></td>
<td><p align=\"center\">&LT;=0</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>q_in.q</h4></p></td>
<td><p align=\"center\">=solutionFlow*q_in.conc</p></td>
<td><p align=\"center\">=-q_out.q</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>q_out.q</h4></p></td>
<td><p align=\"center\">=-q_in.q</p></td>
<td><p align=\"center\">=solutionFlow*q_out.conc</p></td>
</tr>
</table>
<br/>
</html>"));
    end Stream;

    model SubstancePump "Prescribed sunstance molar flow"
      extends Interfaces.OnePortParallel;
      extends Interfaces.ConditionalSubstanceFlow;

    equation
      port_a.q = q;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
                100,100}}), graphics={
            Rectangle(
              extent={{-100,-50},{100,50}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              rotation=360),
            Polygon(
              points={{-80,25},{80,0},{-80,-25},{-80,25}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid,
              rotation=360),
            Text(
              extent={{-150,-20},{150,20}},
              lineColor={0,0,255},
              origin={-10,-76},
              rotation=360,
              textString="%name")}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstancePump;
  end Components;

  package Sensors
    extends Modelica.Icons.SensorsPackage;

    model MolarFlowSensor "Measure of molar flow"

      extends Modelica.Icons.RotationalSensor;
      extends Interfaces.OnePortSerial;

      Modelica.Blocks.Interfaces.RealOutput molarFlowRate(final unit="mol/s") annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-100})));

    equation
      molarFlowRate = port_a.q;

      port_a.u = port_b.u;

     annotation (
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),     Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={
            Line(
              points={{70,-10},{90,-10}},
              color={127,0,127},
              smooth=Smooth.None),
            Line(
              points={{70,10},{90,10}},
              color={127,0,127},
              smooth=Smooth.None),
            Line(
              points={{-90,10},{-70,10}},
              color={127,0,127},
              smooth=Smooth.None),
            Line(
              points={{-90,-10},{-70,-10}},
              color={127,0,127},
              smooth=Smooth.None),
            Text(
              extent={{-31,-5},{28,-64}},
              lineColor={0,0,0},
              textString="dn")}));
    end MolarFlowSensor;

    model MoleFractionSensor "Measure of mole fraction"
      extends Modelica.Icons.RotationalSensor;
      extends Interfaces.PartialSubstanceNoStorage;

      Interfaces.SolutionPort solution annotation (Placement(transformation(
              extent={{-70,-60},{-50,-40}}), iconTransformation(extent={{-70,-60},{-50,
                -40}})));

      Interfaces.SubstanceUsePort port_a "For measure only" annotation (
          Placement(transformation(extent={{-10,-12},{10,8}}),
            iconTransformation(extent={{-10,-12},{10,8}})));
       Modelica.Blocks.Interfaces.RealOutput moleFraction(final unit="1")
        "Mole fraction of the substance"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,0})));

    equation
      port_a.q = 0;

      moleFraction = x;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
                100,100}}), graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
              textString="x"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MoleFractionSensor;

    model MolalitySensor "Measure of molality of the substance"
      extends Modelica.Icons.RotationalSensor;
      extends Interfaces.PartialSubstanceNoStorage;

      parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionPer1kgOfSolvent = 55.508
        "Amount of all particles in the solution per one kilogram of solvent";

      Interfaces.SubstanceUsePort port_a "For measure only" annotation (
          Placement(transformation(extent={{-10,-12},{10,8}}),
            iconTransformation(extent={{-10,-12},{10,8}})));
       Modelica.Blocks.Interfaces.RealOutput molality(final unit="mol/kg")
        "Molality of the substance (amount of substance per mass of solvent)"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,0})));

    protected
      constant Modelica.SIunits.Mass KG=1;
    equation
      port_a.q = 0;

      x=molality*KG / AmountOfSolutionPer1kgOfSolvent;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                            graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
              textString="x"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MolalitySensor;

    model MolarConcentrationSensor "Measure of molarity of the substance"
      extends Modelica.Icons.RotationalSensor;
      extends Interfaces.PartialSubstanceNoStorage;

    parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionInOneLiter = 55.508
        "Amount of all particles in one liter of the solution";

      Interfaces.SubstanceUsePort port_a "For measure only" annotation (
          Placement(transformation(extent={{-10,-12},{10,8}}),
            iconTransformation(extent={{-10,-12},{10,8}})));
       Modelica.Blocks.Interfaces.RealOutput molarConcentration(final unit="mol/m3", displayUnit="mol/l")
        "Molarity of the substance (amount of substance in one liter of whole solution)"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,0})));

    protected
      constant Modelica.SIunits.Volume L=0.001;
    equation
      port_a.q = 0;

      x=molarConcentration*L / AmountOfSolutionInOneLiter;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                            graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
              textString="x"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MolarConcentrationSensor;

    model MassFractionSensor "Measure of mass fraction of the substance"
      extends Modelica.Icons.RotationalSensor;
      extends Interfaces.PartialSubstanceNoStorage;

    parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionInOneKilogram = 55.508
        "Amount of all particles in one kilogram of the solution";

      Interfaces.SubstanceUsePort port_a "For measure only" annotation (
          Placement(transformation(extent={{-10,-12},{10,8}}),
            iconTransformation(extent={{-10,-12},{10,8}})));
       Modelica.Blocks.Interfaces.RealOutput massFraction(final unit="kg/kg")
        "Mass fraction of the substance (mass of the substance per mass of the whole solution)"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,0})));

    equation
      port_a.q = 0;

      x=(massFraction/molarMass) / AmountOfSolutionInOneKilogram;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                            graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
              textString="x"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MassFractionSensor;
  end Sensors;

  package Sources
    extends Modelica.Icons.SourcesPackage;

    model AirSubstance "Substance with defined partial pressure"
      extends Interfaces.PartialSubstance;

      parameter Boolean usePartialPressureInput = false
        "=true, if fixed partial pressure is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Pressure PartialPressure=0
        "Fixed partial pressure if usePartialPressureInput=false"
        annotation (HideResult=true, Dialog(enable=not usePartialPressureInput));

      parameter Modelica.SIunits.Pressure TotalPressure=101325
        "Total pressure of the whole gaseous solution";

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput partialPressure(start=
            PartialPressure, final unit="Pa")=p if usePartialPressureInput
        "Partial pressure of gas = total pressure * gas fraction"
        annotation (HideResult=true,Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.Pressure p "Current partial pressure";

      parameter Modelica.SIunits.Volume Volume = 0.001
        "Volume of gaseous solution";

    equation
      if not usePartialPressureInput then
        p=PartialPressure;
      end if;

      //mole fraction
      x = p / TotalPressure;

      //the solution
      temperature = Temperature;
      pressure = TotalPressure;
      electricPotential = ElectricPotential;
      amountOfSolution = TotalPressure*Volume/(Modelica.Constants.R*Temperature);
      moleFractionBasedIonicStrength = 0;

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={170,255,255},
            fillPattern=FillPattern.Backward),
            Polygon(
              points={{-100,100},{100,-100},{100,100},{-100,100}},
              smooth=Smooth.None,
              fillColor={159,159,223},
              fillPattern=FillPattern.Backward,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Text(
              extent={{0,0},{-100,-100}},
              lineColor={0,0,0},
              textString="P,T"),
            Line(
              points={{-62,0},{56,0}},
              color={191,0,0},
              thickness=0.5),
            Polygon(
              points={{38,-20},{38,20},{78,0},{38,-20}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-100,-102},{104,-126}},
              lineColor={0,0,0},
              textString="%T K")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end AirSubstance;

    model PureSubstance "Constant source of pure substance"
      extends Interfaces.PartialSubstance;

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
        "Electric potential";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
        "Ionic strength";

      parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution = 1
        "Amount of solution";

    equation
      x = 1;

      //the solution
      temperature = Temperature;
      pressure = Pressure;
      electricPotential = ElectricPotential;
      moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;
      amountOfSolution = AmountOfSolution;

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillColor={107,45,134},
              fillPattern=FillPattern.Backward),
            Text(
              extent={{10,8},{-90,-92}},
              lineColor={0,0,0},
              textString="pure"),
            Line(
              points={{-62,0},{56,0}},
              color={191,0,0},
              thickness=0.5),
            Polygon(
              points={{38,-20},{38,20},{78,0},{38,-20}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-104,-76},{100,-100}},
              lineColor={0,0,0},
              textString="%T K")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PureSubstance;

    model PureElectricParticle
      "Constant source of pure particle driven by electric port"
      extends Interfaces.PartialSubstance;

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
        "Ionic strength";

      Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
          Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
              extent={{-110,-10},{-90,10}})));

      parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution = 1
        "Amount of solution";

    equation
      //electric
      pin.v = electricPotential;
      pin.i = -Modelica.Constants.F*port_a.q;

      //pure substance
      x = 1;

      //the solution
      temperature = Temperature;
      pressure = Pressure;
      moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;
      amountOfSolution = AmountOfSolution;

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillColor={107,45,134},
              fillPattern=FillPattern.Backward),
            Text(
              extent={{10,8},{-90,-92}},
              lineColor={0,0,0},
              textString="pure"),
            Line(
              points={{-62,0},{56,0}},
              color={191,0,0},
              thickness=0.5),
            Polygon(
              points={{38,-20},{38,20},{78,0},{38,-20}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-104,-76},{100,-100}},
              lineColor={0,0,0},
              textString="%T K")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PureElectricParticle;

    model AmbientMolality "Constant source of substance molality"
      extends Interfaces.PartialSubstance;

       parameter Real Molality(final unit="mol/kg") = 1e-8
        "Fixed molality of the substance if useMolalityInput=false"
        annotation (HideResult=true, Dialog(enable=not useMolalityInput));

       parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution = 55.508
        "Amount of all particles in the solution per one kilogram of solvent";

        parameter Boolean useMolalityInput = false
        "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput molalityInput(start=Molality,final unit="mol/kg")=n/KG if
           useMolalityInput
        annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.AmountOfSubstance n "Current amount of the substance";

    protected
      constant Modelica.SIunits.Mass KG=1;
    equation
       if not useMolalityInput then
         n=Molality*KG;
       end if;

      x = n/AmountOfSolution;

      //solution properties at the port
      temperature = Temperature;
      pressure = Pressure;
      electricPotential = ElectricPotential;
      moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;
      amountOfSolution = AmountOfSolution;

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillColor={107,45,134},
              fillPattern=FillPattern.Backward),
            Line(
              points={{-62,0},{56,0}},
              color={191,0,0},
              thickness=0.5),
            Polygon(
              points={{38,-20},{38,20},{78,0},{38,-20}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-104,-76},{100,-100}},
              lineColor={0,0,0},
              textString="%T K"),
            Text(
              extent={{94,-4},{-94,-78}},
              lineColor={0,0,0},
              textString="molality")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end AmbientMolality;

    model AmbientConcentration "Constant source of molar concentration"
       extends Interfaces.PartialSubstance;

       parameter Real MolarConcentration(final unit="mol/m3", displayUnit="mol/l") = 1e-8
        "Fixed molarity of the substance if useMolarityInput=false"
        annotation (HideResult=true, Dialog(enable=not useMolarityInput));

       parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution = 55.508
        "Amount of all particles in the solution one liter of solvent";

        parameter Boolean useMolarityInput = false
        "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

       parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput molarConcentrationInput(start=MolarConcentration,final unit="mol/m3", displayUnit="mol/l")=n/L if
           useMolarityInput
        annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.AmountOfSubstance n "Current amount of the substance";

    protected
      constant Modelica.SIunits.Volume L=0.001;
    equation
       if not useMolarityInput then
         n=MolarConcentration*L;
       end if;

      x = n/AmountOfSolution;

      //solution properties at the port
      temperature = Temperature;
      pressure = Pressure;
      electricPotential = ElectricPotential;
      moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;
      amountOfSolution = AmountOfSolution;

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillColor={107,45,134},
              fillPattern=FillPattern.Backward),
            Text(
              extent={{94,92},{-94,18}},
              lineColor={0,0,0},
              textString="molarity"),
            Line(
              points={{-62,0},{56,0}},
              color={191,0,0},
              thickness=0.5),
            Polygon(
              points={{38,-20},{38,20},{78,0},{38,-20}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-104,-76},{100,-100}},
              lineColor={0,0,0},
              textString="%T K")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end AmbientConcentration;

    model AmbientMoleFraction "Constant source of substance mole fraction"
         extends Interfaces.PartialSubstance;

       parameter Modelica.SIunits.MoleFraction MoleFraction = 1e-8
        "Fixed mole fraction of the substance if useMoleFractionInput=false"
        annotation (HideResult=true, Dialog(enable=not useMoleFractionInput));

       parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution = 55.508
        "Amount of all reacting particles in the solution";

        parameter Boolean useMoleFractionInput = false
        "Is mole fraction of the substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput moleFractionInput(
        final unit="mol/mol",
        start=MoleFraction)=x if
           useMoleFractionInput annotation (HideResult=true, Placement(transformation(
              extent={{-120,-20},{-80,20}})));

    equation
       if not useMoleFractionInput then
         x=MoleFraction;
       end if;

      //solution properties at the port
      temperature = Temperature;
      pressure = Pressure;
      electricPotential = ElectricPotential;
      moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;
      amountOfSolution = AmountOfSolution;

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillColor={107,45,134},
              fillPattern=FillPattern.Backward),
            Line(
              points={{-62,0},{56,0}},
              color={191,0,0},
              thickness=0.5),
            Polygon(
              points={{38,-20},{38,20},{78,0},{38,-20}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-104,-76},{100,-100}},
              lineColor={0,0,0},
              textString="%T K"),
            Text(
              extent={{94,-4},{-94,-78}},
              lineColor={0,0,0},
              textString="n")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end AmbientMoleFraction;

    model SubstanceInflow "Molar pump of substance to system"
      extends Interfaces.ConditionalSubstanceFlow;

      Interfaces.SubstanceUsePort port_b "Outflow"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    equation
      port_b.q = - q;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
                100,100}}), graphics={
            Rectangle(
              extent={{-100,-42},{100,40}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-48,20},{50,0},{-48,-21},{-48,20}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-82,-82},{90,-58}},
              textString="%name",
              lineColor={0,0,255})}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceInflow;

    model SubstanceOutflow "Molar pump of substance out of system"
      extends Interfaces.ConditionalSubstanceFlow;

      Interfaces.SubstanceUsePort port_a "Inflow"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

    equation
      port_a.q = q;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
                100,100}}), graphics={
            Rectangle(
              extent={{-100,-42},{100,40}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-48,20},{50,0},{-48,-21},{-48,20}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-82,-82},{90,-58}},
              textString="%name",
              lineColor={0,0,255})}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceOutflow;

    model Clearance "Physiological Clearance"
     extends Interfaces.ConditionalSolutionFlow(final SolutionFlow=Clearance/K);
     extends Interfaces.PartialSubstanceNoStorage;

      parameter Modelica.SIunits.VolumeFlowRate Clearance=0
        "Physiological clearance of the substance if useSolutionFlowInput=false"
        annotation (HideResult=true, Dialog(enable=not useSolutionFlowInput));

      parameter Real K(unit="1")=1
        "Coefficient such that Clearance = K*solutionFlow";

      Modelica.SIunits.MolarFlowRate molarClearance "Current molar clearance";

    equation
      molarClearance = q*K;

      port_a.q = molarClearance * x;

      assert(molarClearance>=-Modelica.Constants.eps, "Clearance can not be negative!");

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                            graphics={
            Rectangle(
              extent={{-100,-100},{100,50}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{80,25},{-80,0},{80,-25},{80,25}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,-90},{150,-50}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-100,-30},{100,-50}},
              lineColor={0,0,0},
              textString="K=%K")}),        Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end Clearance;

    model Degradation "Degradation of the substance"
      extends Interfaces.PartialSubstanceNoStorage;

      parameter Physiolibrary.Types.Time HalfTime
        "Degradation half time. The time after which will remain half of initial concentration in the defined volume when no other generation, clearence and degradation exist.";

    equation
      port_a.q = (Modelica.Math.log(2)/HalfTime)*x*amountOfSolution;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                            graphics={
            Rectangle(
              extent={{-100,-100},{100,58}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{64,26},{-78,0},{64,-26},{64,26}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-148,-82},{152,-42}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-100,54},{100,28}},
              lineColor={0,0,0},
              textString="t1/2 = %HalfTime s"),
            Polygon(
              points={{54,24},{54,-24},{44,-22},{44,22},{54,24}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{30,20},{30,-20},{20,-18},{20,18},{30,20}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{8,16},{8,-16},{-2,-14},{-2,14},{8,16}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-12,12},{-12,-12},{-22,-10},{-22,10},{-12,12}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-34,8},{-34,-8},{-44,-6},{-44,6},{-34,8}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-56,4},{-56,-4},{-66,-2},{-66,2},{-56,4}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid)}),
        Documentation(revisions="<html>
<table>
<tr>
<td>Author:</td>
<td>Marek Matejak</td>
</tr>
<tr>
<td>Copyright:</td>
<td>In public domains</td>
</tr>
<tr>
<td>By:</td>
<td>Charles University, Prague</td>
</tr>
<tr>
<td>Date of:</td>
<td>2013-2015</td>
</tr>
</table>
</html>"));
    end Degradation;

    model Buffer
      "Source of substance to reach linear dependence between concentration and electrochemical potential"
         extends Interfaces.PartialSubstance;

       parameter Modelica.SIunits.MoleFraction xBuffered_start=1e-7
        "Initial value of mole fraction of the buffered substance";

       parameter Real BufferValue(final unit="1") = 1
        "Fixed buffer value (slope between x and -log10(x)) if useBufferValueInput=false"
        annotation (HideResult=true, Dialog(enable=not useMoleFractionInput));

       parameter Boolean useBufferValueInput = false
        "Is buffer value of the substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

        Real bufferValue(final unit="1");

      Modelica.Blocks.Interfaces.RealInput bufferValueInput(
        final unit="mol/mol",
        start=BufferValue)=bufferValue if
           useBufferValueInput annotation (HideResult=true, Placement(transformation(
              extent={{-120,-20},{-80,20}})));

       parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution = 55.508
        "Amount of all particles in the buffering solution";

      parameter Modelica.SIunits.Time Tau = 1
        "Time constant for other scaling of the buffering rate";

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
        "Electric potential";

      Modelica.SIunits.AmountOfSubstance nBuffered;
      Modelica.SIunits.MoleFraction xBuffered;
    initial equation
      xBuffered = xBuffered_start;
    equation
      if not useBufferValueInput then
        bufferValue = BufferValue;
      end if;

      der(nBuffered) = port_a.q;
      xBuffered = nBuffered/AmountOfSolution;
      port_a.q = (AmountOfSolution/Tau)*(-xBuffered -log10(x)*bufferValue);

      //solution properties at the port
      temperature = Temperature;
      pressure = Pressure;
      electricPotential = ElectricPotential;
      moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;
      amountOfSolution = AmountOfSolution;

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.CrossDiag),
            Line(
              points={{-62,0},{56,0}},
              color={191,0,0},
              thickness=0.5),
            Polygon(
              points={{38,-20},{38,20},{78,0},{38,-20}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-104,-76},{100,-100}},
              lineColor={0,0,0},
              textString="%T K"),
            Text(
              extent={{94,-4},{-94,-78}},
              lineColor={0,0,0},
              textString="n")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Buffer;
  end Sources;

  package Interfaces
    extends Modelica.Icons.InterfacesPackage;

    package SubstanceModel "Base substance model"

       record SubstanceData "Base substance data"

       parameter Modelica.SIunits.MolarMass MolarWeight(displayUnit="kDa")=0
          "Molar weight of the substance in kg/mol or kDa";

       parameter Modelica.SIunits.ChargeNumberOfIon z=0
          "Charge number of the substance (e.g. 0..uncharged, -1..electron, +2..Ca^2+)";

       parameter Modelica.SIunits.MolarEnergy DfH(displayUnit="kJ/mol")=0
          "Enthalpy of formation of the substance in the selected state";
       parameter Modelica.SIunits.MolarEnergy DfG_25degC(displayUnit="kJ/mol")=0
          "Gibbs enerfy of formation at 25°C of the substance in the selected state";

       parameter String References[:]={""}
          "References of these thermodynamical values";

        annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
       end SubstanceData;

     replaceable function activityCoefficient
        "Return activity coefficient of the substance in the solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Real gamma "Activity Coefficient";
     algorithm
         gamma := 1;
     end activityCoefficient;

     replaceable function chargeNumberOfIon
        "Return charge number of the substance in the solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.ChargeNumberOfIon z "Charge number of ion";
     algorithm
        z := substanceData.z;
     end chargeNumberOfIon;

     replaceable function molarEnthalpy "Molar enthalpy of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarEnthalpy molarEnthalpy "Molar enthalpy";
     algorithm
         molarEnthalpy := substanceData.DfH + Modelica.Constants.F*substanceData.z*v;
     end molarEnthalpy;

     replaceable function molarEntropy "Molar entropy of the substance"
        extends Modelica.Icons.Function;
        input Modelica.SIunits.ChemicalPotential u
          "Electro-chemical potential of the substance";
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarEntropy molarEntropy "Molar entropy";
     algorithm
         molarEntropy :=  (u - molarEnthalpy(substanceData,T,p,v,I))/T;
     end molarEntropy;

     replaceable function u0
        "Chemical part of electro-chemical potential of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.ChemicalPotential u0 "Base chemical potential";
     algorithm
         u0 := substanceData.DfH - T*((substanceData.DfH-substanceData.DfG_25degC)/298.15);
     end u0;

     replaceable function uPure
        "Electro-chemical potential of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.ChemicalPotential uPure
          "Base electro-chemical potential";
     algorithm
         uPure := u0(substanceData,T,p,v,I) + Modelica.Constants.F*substanceData.z*v;
     end uPure;

     replaceable function molarMass "Molar mass of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarMass molarMass "Molar mass";
     algorithm
         molarMass := substanceData.MolarWeight; //ideal gas
         //incompressible: molarVolume := constant;
     end molarMass;

     //not needed for isobaric processes:
     replaceable function molarVolume "Molar volume of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
          "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarVolume molarVolume "Molar volume";
     algorithm
         molarVolume := Modelica.Constants.R*T/p; //ideal gas
         //incompressible: molarVolume := constant;
     end molarVolume;

      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceModel;

    connector ChemicalPort
      "Electro-chemical potential and molar change of the substance in the solution"

      Modelica.SIunits.ChemicalPotential u
        "Electro-chemical potential of the substance in the solution";

      flow Modelica.SIunits.MolarFlowRate q "Molar change of the substance";

      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Definition of electro-chemical potential of the substance:</p>
<p><b>u<sub>e-ch</sub>(x,T,V) = u&deg;(T) + R*T*ln(gamma*x) + z*F*V</b></p>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>V .. eletric potential of the substance</p>
<p>z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
<p>R .. gas constant</p>
<p>F .. Faraday constant</p>
<p>gamma .. activity coefficient</p>
<p>u&deg;(T) .. chemical potential of pure substance</p>
<p>DfG(T) .. free Gibbs energy of formation of the substance at current temperature T. </p>
<p>DfH .. free enthalpy of formation of the substance</p>
<p>DfS .. free entropy of formation of the substance </p>
<p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
</html>"));
    end ChemicalPort;

    connector SubstanceDefinitionPort
      "Electro-chemical potential and molar flow of the substance in the solution"
      extends ChemicalPort;

    /*  
 //substance properties (expressed from sunstance definition and current state of the solution)
  output Modelica.SIunits.MolarMass molarWeight "Molar weight of the substance";

  output Modelica.SIunits.ChargeNumberOfIon z "Charge number of the substance";

  //substance properties dependent on solution
  output Modelica.SIunits.MoleFraction x 
    "Mole fraction of the substance in the solution";

  output Modelica.SIunits.MolarEnthalpy molarEnthalpy(displayUnit="kJ/mol") 
    "Molar enthalpy of the substance";

  output Modelica.SIunits.MolarEntropy molarEntropy 
    "Molar entropy of the substance";

  output Modelica.SIunits.MolarVolume molarVolume 
    "Molar volume of the substance";

  output Modelica.SIunits.MolarEnergy u0 
    "Chemical potential of the pure substance";

  output Modelica.SIunits.MolarEnergy uPure 
    "Electro-Chemical potential of the pure substance";

  output Modelica.SIunits.ActivityCoefficient activityCoefficient 
    "Activity coefficient of the substance";

  //solution properties
  output Modelica.SIunits.Temperature temperature(displayUnit="degC") 
    "Temperature of the solution";

  output Modelica.SIunits.Pressure pressure(displayUnit="bar") 
    "Pressure of the solution";

  output Modelica.SIunits.ElectricPotential electricPotential 
    "Total electric potential of the solution";

  output Modelica.SIunits.AmountOfSubstance amountOfSolution 
    "Total amount of all particles in the solution";
    */

    annotation (
        defaultComponentName="port_a",
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                100}}),     graphics={Rectangle(
              extent={{-20,10},{20,-10}},
              lineColor={158,66,200},
              lineThickness=1),       Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={158,66,200},
              fillColor={158,66,200},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}),
            graphics={Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={158,66,200},
              fillColor={158,66,200},
              fillPattern=FillPattern.Solid,
              lineThickness=1),
       Text(extent = {{-160,110},{40,50}}, lineColor={172,72,218},   textString = "%name")}),
        Documentation(info="<html>
<p>Chemical port with internal definition of the substance inside the component. </p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceDefinitionPort;

    connector SubstanceUsePort
      "Electro-chemical potential and molar flow of the substance in the solution"
      extends ChemicalPort;

    /*
  //substance properties (read only; user does not need to set them again)
  input Modelica.SIunits.MolarMass molarWeight "Molar weight of the substance";

  input Modelica.SIunits.ChargeNumberOfIon z "Charge number of the substance";

  //substance properties dependent on solution
  input Modelica.SIunits.MoleFraction x 
    "Mole fraction of the substance in the solution";

  input Modelica.SIunits.MolarEnthalpy molarEnthalpy(displayUnit="kJ/mol") 
    "Molar enthalpy of the substance in the solution";

  input Modelica.SIunits.MolarEntropy molarEntropy 
    "Molar entropy of the substance in the solution";

  input Modelica.SIunits.MolarVolume molarVolume 
    "Molar volume of the substance";

  input Modelica.SIunits.MolarEnergy u0 
    "Chemical potential of the pure substance";

  input Modelica.SIunits.MolarEnergy uPure 
    "Electro-Chemical potential of the pure substance";

  input Modelica.SIunits.ActivityCoefficient activityCoefficient 
    "Activity coefficient of the substance in the solution";

  //solution properties
  input Modelica.SIunits.Temperature temperature(displayUnit="degC") 
    "Temperature of the solution";

  input Modelica.SIunits.Pressure pressure(displayUnit="bar") 
    "Pressure of the solution";

  input Modelica.SIunits.ElectricPotential electricPotential 
    "Total electric potential of the solution";

  input Modelica.SIunits.AmountOfSubstance amountOfSolution 
  "Total amount of all particles in the solution";
  */

    annotation (
        defaultComponentName="port_a",
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                100}}),     graphics={Rectangle(
              extent={{-20,10},{20,-10}},
              lineColor={158,66,200},
              lineThickness=1),       Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={158,66,200},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}),
            graphics={Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={158,66,200},
              fillColor={158,66,200},
              fillPattern=FillPattern.Solid,
              lineThickness=1),
       Text(extent = {{-160,110},{40,50}}, lineColor={172,72,218},   textString = "%name")}),
        Documentation(info="<html>
<p>Chemical port with external definition of the substance outside the component.</p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceUsePort;

    partial model PartialSubstance

      Interfaces.SubstanceDefinitionPort port_a
        "The substance with prescribed partial pressure"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      parameter substanceModel.SubstanceData substanceData
         annotation (choicesAllMatching = true);

      replaceable package substanceModel = Interfaces.SubstanceModel   constrainedby
        Interfaces.SubstanceModel
        "Substance model: Molar Weight, Enthalpy, Gibbs energy,... from substance data"
         annotation (choicesAllMatching = true);

      Modelica.SIunits.MoleFraction x "Mole fraction of the substance";

    protected
      Modelica.SIunits.ActivityCoefficient gamma
        "Activity coefficient of the substance";

      Modelica.SIunits.ChargeNumberOfIon z "Charge number of ion";

      Modelica.SIunits.Temperature temperature(start=298.15)
        "Temperature of the solution";

      Modelica.SIunits.Pressure pressure(start=101325)
        "Pressure of the solution";

      Modelica.SIunits.ElectricPotential electricPotential(start=0)
        "Electric potential of the solution";

      Modelica.SIunits.MoleFraction moleFractionBasedIonicStrength(start=0)
        "Ionic strength of the solution";

      Modelica.SIunits.AmountOfSubstance amountOfSolution
        "Amount of all solution particles";

      Modelica.SIunits.MolarMass molarMass "Molar mass of the substance";

      Modelica.SIunits.MolarEnthalpy molarEnthalpy
        "Molar enthalpy of the substance";

      Modelica.SIunits.MolarEntropy molarEntropy
        "Molar entropy of the substance";

      Modelica.SIunits.ChemicalPotential u0
        "Chemical potential of the pure substance";

      Modelica.SIunits.ChemicalPotential uPure
        "Electro-Chemical potential of the pure substance";

      Modelica.SIunits.MolarVolume molarVolume "Molar volume of the substance";

    equation
      //define the substance
     /* port_a.x = x;
  port_a.activityCoefficient = substanceModel.activityCoefficient(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
  port_a.molarWeight = substanceModel.molarMass(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
  port_a.z =  substanceModel.chargeNumberOfIon(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);

  port_a.molarEnthalpy = substanceModel.molarEnthalpy(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
  port_a.molarEntropy = substanceModel.molarEntropy(port_a.u,substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
  port_a.u0 = substanceModel.u0(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
  port_a.uPure = substanceModel.uPure(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
  port_a.molarVolume = substanceModel.molarVolume(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);

  //the solution
  port_a.temperature = temperature;
  port_a.pressure = pressure;
  port_a.electricPotential = electricPotential;
  port_a.amountOfSolution = amountOfSolution;
  */

      gamma = substanceModel.activityCoefficient(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
      z = substanceModel.chargeNumberOfIon(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
      molarMass = substanceModel.molarMass(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);

      molarEnthalpy = substanceModel.molarEnthalpy(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
      molarEntropy = substanceModel.molarEntropy(port_a.u,substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
      u0 = substanceModel.u0(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
      uPure = substanceModel.uPure(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
      molarVolume = substanceModel.molarVolume(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);

      //electro-chemical potential of the substance in the solution
      port_a.u = substanceModel.u0(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength)
       + Modelica.Constants.R*temperature*log(gamma*x)
       + z*Modelica.Constants.F*electricPotential;

      annotation (
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PartialSubstance;

    partial model OnePortParallel
      "Partial molar flow beween two substance definitions"

      SubstanceUsePort port_a annotation (Placement(transformation(extent={{-110,-10},
                {-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));
      SubstanceUsePort port_b annotation (Placement(transformation(extent={{90,-10},
                {110,10}}), iconTransformation(extent={{90,-10},{110,10}})));
    equation
      port_a.q + port_b.q = 0;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end OnePortParallel;

    partial model OnePortSerial
      "Partial transfer of substance from substance definition component to another transfer component (such as MolarFlowSensor)"

      SubstanceUsePort port_a annotation (Placement(transformation(extent={{-110,
                -10},{-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));
      SubstanceDefinitionPort port_b
        annotation (Placement(transformation(extent={{90,-10},{110,10}}),
            iconTransformation(extent={{90,-10},{110,10}})));
    equation
      port_a.q + port_b.q = 0;
    /*
  port_b.molarWeight = port_a.molarWeight;
  port_b.z = port_a.z;
  port_b.x = port_a.x;
  port_b.molarEnthalpy = port_a.molarEnthalpy;
  port_b.molarEntropy = port_a.molarEntropy;
  port_b.molarVolume = port_a.molarVolume;
  port_b.u0 = port_a.u0;
  port_b.uPure = port_a.uPure;
  port_b.activityCoefficient = port_a.activityCoefficient;
  port_b.temperature = port_a.temperature;
  port_b.pressure = port_a.pressure;
  port_b.electricPotential = port_a.electricPotential;
  port_b.amountOfSolution = port_a.amountOfSolution;
*/
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end OnePortSerial;

    partial model ConditionalSolutionFlow
      "Input of solution molar flow vs. parametric solution molar flow"

      parameter Boolean useSolutionFlowInput = false
        "Is solution flow an input?"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.VolumeFlowRate SolutionFlow=0
        "Volume flow rate of the solution if useSolutionFlowInput=false" annotation (
          HideResult=true, Dialog(enable=not useSolutionFlowInput));

      parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L=55.508
        "The amount of all particles in one liter of the solution";

      Modelica.Blocks.Interfaces.RealInput solutionFlow(start=SolutionFlow, final unit="m3/s")=
         q*OneLiter/AmountOfSolutionIn1L if useSolutionFlowInput
         annotation ( HideResult=true, Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,40}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,70})));

      Modelica.SIunits.MolarFlowRate q "Current molar solution flow";

    protected
     constant Modelica.SIunits.Volume OneLiter=0.001 "One liter";

    equation
      if not useSolutionFlowInput then
        q*OneLiter/AmountOfSolutionIn1L = SolutionFlow;
      end if;

    end ConditionalSolutionFlow;

    partial model ConditionalSubstanceFlow
      "Input of substance molar flow vs. parametric substance molar flow"

      parameter Boolean useSubstanceFlowInput = false
        "Is substance flow an input?"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.MolarFlowRate SubstanceFlow=0
        "Volumetric flow of Substance if useSubstanceFlowInput=false" annotation (
          HideResult=true, Dialog(enable=not useSubstanceFlowInput));

      Modelica.Blocks.Interfaces.RealInput substanceFlow(start=SubstanceFlow, final unit="mol/s")=q if
           useSubstanceFlowInput
           annotation (HideResult=true,
           Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={40,40})));

      Modelica.SIunits.MolarFlowRate q "Current Substance flow";
    equation
      if not useSubstanceFlowInput then
        q = SubstanceFlow;
      end if;

    end ConditionalSubstanceFlow;

    connector SolutionPort
      "Interation of properties from all substances of the solution"

      //heat port
      Modelica.SIunits.Temperature T "Temperature of the solution";
      flow Modelica.SIunits.EnthalpyFlowRate dH
        "Enthalpy change of the solution";

      //entropy
      Modelica.SIunits.Entropy S "Free entropy";
      flow Modelica.SIunits.EntropyFlowRate dS "Entropy change of the solution";

      //free Gibbs energy
      Modelica.SIunits.Energy G "Free Gibbs energy";
      flow Modelica.SIunits.EnergyFlowRate dG
        "Gibbs energy change of the solution";

      //amount of substances
      Modelica.SIunits.AmountOfSubstance n "Amount of the solution";
      flow Modelica.SIunits.MolarFlowRate dn "Molar change of the solution";

      //electric port
      Modelica.SIunits.ElectricPotential v "Electric potential in the solution";
      flow Modelica.SIunits.ElectricCurrent i "Change of electric charge";

      //ionic strength
      Modelica.SIunits.MoleFraction I
        "Mole fraction based ionic strength of the solution";
      flow Modelica.SIunits.MolarFlowRate dI
        "Change of mole-fraction based ionic strength of the solution";

      //hydraulic port
      Modelica.SIunits.Pressure p "Pressure of the solution";
      flow Modelica.SIunits.VolumeFlowRate dV "Volume change of the solution";

      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Solution port integrates all substances of the solution:</p>
<p>Such as if there are conected together with electric port, thermal port and with port composed with the amont of substance and molar change of substance.</p>
</html>"), Icon(graphics={            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={158,66,200},
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid)}));
    end SolutionPort;

    partial model PartialSolution
      "Chemical solution as homogenous mixture of the substances (only pressure and electric potential are not defined)"

      parameter Modelica.SIunits.MolarHeatCapacity Cp_start = 75.4
        "Initial heat capacity of the solution at constant pressure"
        annotation (Dialog(group="Initialization"));

      parameter Modelica.SIunits.AmountOfSubstance amountOfSolution_start=55.508
        "Initial amount of the solution (default is amount of pure water in 1kg)"
         annotation (Dialog(group="Initialization"));

      parameter Modelica.SIunits.Volume volume_start=1/0.997
        "Initial volume of the solution (default is volume of pure water in 1kg at 25degC)"
         annotation (Dialog(group="Initialization"));

      parameter Modelica.SIunits.Temperature temperature_start=298.15
        "Initial temperature of the solution"
         annotation (Dialog(group="Initialization"));

      parameter Modelica.SIunits.ElectricCharge electricCharge_start=0
        "Initial electric charge of the solution"
         annotation (Dialog(group="Initialization"));

      parameter Modelica.SIunits.MoleFraction ionicStrength_start=0
        "Initial ionic strength (mole fraction based) of the solution"
         annotation (Dialog(group="Initialization"));

      Modelica.SIunits.AmountOfSubstance amountOfSolution(start=amountOfSolution_start)
        "Current amount of all substances in the solution";

      Modelica.SIunits.ElectricCharge charge(start=electricCharge_start)
        "Current amount of all substances in the solution";

      Modelica.SIunits.Enthalpy freeEnthalpy(start=temperature_start*Cp_start*amountOfSolution_start)
        "Free enthalpy of the solution";

      Interfaces.SolutionPort solution "Solution nonflows and flows"
                                      annotation (Placement(
            transformation(extent={{-80,-80},{-60,-60}}),iconTransformation(extent={{-2,-102},{2,-98}})));

      //for debuging only:
      Modelica.SIunits.MolarHeatCapacity Cp=freeEnthalpy/(solution.T*amountOfSolution)
        "Current heat capacity of the solution";

      //Valid only, if you are sure with molarVolume calculation of the substances.
      Modelica.SIunits.Volume volume(start=volume_start)
        "Current volume of the solution (Valid only, if you are sure with molarVolume calculation of the substances)";

      Modelica.SIunits.HeatFlowRate heatFromEnvironment;

    //  Real lnn(stateSelect=StateSelect.prefer) "ln(solution.n)";
    //  Real lncharge(stateSelect=StateSelect.prefer) "ln(charge)";
    //  Real lnIn(stateSelect=StateSelect.prefer) "ln(solution.I*amountOfSolution)";
    //  Real lnvolume(stateSelect=StateSelect.prefer) "ln(volume)";

    initial equation
      amountOfSolution = amountOfSolution_start;
      volume = volume_start;
      charge = electricCharge_start;
      freeEnthalpy = temperature_start*Cp_start*amountOfSolution_start;
      solution.S = Cp_start*amountOfSolution_start;
      solution.G = 0;
      solution.I = ionicStrength_start;
    equation

      //amount of substances
      der(amountOfSolution) = solution.dn;
      //der(lnn) = solution.dn/amountOfSolution;  amountOfSolution = exp(lnn);
      solution.n = amountOfSolution;

      //heat
      der(solution.S) = solution.dS;
      der(freeEnthalpy) = solution.dH + heatFromEnvironment;
      der(solution.G) = solution.dG;
      solution.T = (freeEnthalpy-solution.G)/solution.S;

      //electric
      der(charge) = solution.i;
      //der(lncharge) = solution.i/charge; charge = exp(lncharge);

      //ionic strength (mole fraction based)
      der(solution.I*amountOfSolution) = solution.dI;
      //der(lnIn) = solution.dI/(solution.I*amountOfSolution); (solution.I*amountOfSolution) = exp(lnIn);

      //isobaric
      der(volume) = solution.dV;
      //der(lnvolume) = solution.dV/volume; volume = exp(lnvolume);

                                                                                                          annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false, initialScale=1, extent={{-100,-100},{100,100}}),
            graphics={
            Line(
              points={{-100,100},{-100,72},{-100,72}},
              color={127,0,127},
              smooth=Smooth.None),
            Line(
              points={{100,100},{100,80},{100,70}},
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p>2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<h4>amountOfSubstances = &int; MolarFlows</h4>
<h4>electricCharge = &int; ElectricCurrents</h4>
<h4>freeEnthalpy = &int; EnthalpyChanges</h4>
<h4>freeEntropy = &int; EntropyChanges</h4>
<h4>freeGibbsEnergy = &int; GibbsEnergyChanges</h4>
<h4>electricEnergy = &int; ElectricPowers</h4>
<p>Integration of all substances together into one homogenous mixture - the solution.</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end PartialSolution;

    partial model PartialSubstanceInSolution
      "Substance properties for components, where the substance is connected with the solution"
      extends PartialSubstance;

      SolutionPort            solution
        "To connect substance with solution, where is pressented"                                annotation (Placement(transformation(
              extent={{-70,-110},{-50,-90}}),iconTransformation(extent={{-70,-110},{
                -50,-90}})));
    equation
      //global properties of the solution
      temperature = solution.T;
      pressure = solution.p;
      electricPotential = solution.v;
      amountOfSolution = solution.n;
      moleFractionBasedIonicStrength = solution.I;
    end PartialSubstanceInSolution;

    partial model PartialSubstanceNoStorage
      "Substance properties for components, where the substance is not accumulated"
      extends PartialSubstanceInSolution;

    equation
      //solution is not changed by the non-accumulating components
      //duality allows to change it only by accumulation components
      solution.dH = 0;
      solution.dS = 0;
      solution.dG = 0;
      solution.dn = 0;
      solution.i = 0;
      solution.dI = 0;
      solution.dV = 0;
    end PartialSubstanceNoStorage;
  end Interfaces;

  package Icons "Icons for chemical models"
    //extends Modelica.Icons.IconsPackage;
    extends Modelica.Icons.Package;

    partial class Diffusion

      annotation (Icon(graphics={Bitmap(extent={{-100,100},{100,-100}}, fileName=
                  "modelica://Chemical/Resources/Icons/diffusion.png")}));

    end Diffusion;

    class Substance

        annotation ( Icon(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics={Bitmap(extent={{-100,100},{100,-100}}, fileName=
                  "modelica://Chemical/Resources/Icons/Concentration.png")}));
    end Substance;

    class Speciation

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics={Bitmap(extent={{-100,100},{100,-100}}, fileName=
                  "modelica://Chemical/Resources/Icons/Speciation.png")}));
    end Speciation;

    class GasSolubility

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Bitmap(extent={{-100,100},{100,-100}},
                fileName=
                  "modelica://Chemical/Resources/Icons/GasSolubility.png")}));
    end GasSolubility;

    class Membrane

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Bitmap(extent={{-100,100},{100,-100}},
                fileName="modelica://Chemical/Resources/Icons/membrane.png")}));
    end Membrane;

    class EnzymeKinetics

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Bitmap(extent={{-80,84},{86,-26}},
                fileName=
                  "modelica://Chemical/Resources/Icons/MichaelisMenten.png")}));
    end EnzymeKinetics;

    annotation (Documentation(revisions=""));
  end Icons;

  annotation (
preferredView="info",
version="1.0.0-alpha",
versionBuild=1,
versionDate="2015-04-24",
dateModified = "2015-04-24 17:14:41Z",
uses(Modelica(version="3.2.1"), Physiolibrary(version="2.3.0-beta")),
  Documentation(revisions="<html>
<p>Licensed by Marek Matejak under the Modelica License 2</p>
<p>Copyright &copy; 2008-2015, Marek Matejak, Charles University in Prague.</p>
<p><br><i>This Modelica package is&nbsp;<u>free</u>&nbsp;software and the use is completely at&nbsp;<u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see&nbsp;<a href=\"modelica://Physiolibrary.UsersGuide.ModelicaLicense2\">UsersGuide.ModelicaLicense2</a>&nbsp;or visit&nbsp;<a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>", info="<html>
<p>In physiology books, chapters about chemical substances are organized by their types. The main reason for this is that each substance in the human body is regulated in a different way. For example the regulation of sodium is different from the regulation of potassium, and from the regulation of glucose, and so on. This view leads to the idea of having separate models of each substance. The origin of different flows and regulations is the (cellular) membrane. Water and solutions can cross it in different directions at the same time. Crossings occur for different reasons: water is driven mostly by osmotic gradients, electrolytes are driven by charge to reach Donnan&apos;s equilibrium, and some solutes can even be actively transported against their concentration or electrical gradients. And all this is specifically driven from the higher levels by neural and hormonal responses.&nbsp; </p>
<p>In Physiolibrary flows and fluxes of solutes are supported mostly by the Chemical package. All parts inside this Chemical package use the connector ChemicalPort, which defines the molar concentration and molar flow/flux rate of one solute. This is the supporting infrastructure for modeling membrane diffusion, accumulations of substances, reversal chemical reactions, Henry&apos;s law of gas solubility, dilution with additional solvent flow, membrane reabsorption, chemical degradation and physiological clearance. </p>
</html>"));

end Chemical;
