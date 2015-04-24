within ;
package Chemical
  "Chemical library (reactions, diffusions, semipermeable membranes, gas dissolutions, electrochemical cells, ...)"
 extends Modelica.Icons.Package;

  package Examples
    "Examples that demonstrate usage of the Pressure flow components"
  extends Modelica.Icons.ExamplesPackage;

    package Substances "Definitions of substances"
        extends Modelica.Icons.Package;

      package Silver_solid =
          Interfaces.BaseSubstance (
            MolarWeight=0.1078682,
            z=0,
            DfH=0,
            DfG_25degC=0,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"}) "Ag(s)";
      package Silver_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.1078682,
            z=1,
            DfH=105900,
            DfG_25degC=77100,
            References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "Ag+(aq)";

      package SilverChloride_solid =
          Interfaces.BaseSubstance (
            MolarWeight=0.14332,
            z=0,
            DfH=-127030,
            DfG_25degC=-109720,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "AgCl(s)";

       package Calcium_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.0401,
          z=2,
          DfH=-542960,
            DfG_25degC=-542960-298.15*(33.67),
            References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
        "}) "Ca++(aq)";

      package Chloride_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.03545,
          z=-1,
          DfH=-167460,
            DfG_25degC=-131170,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "Cl-(aq)";

      package CarbonDioxide_gas =
          Interfaces.BaseSubstance (
            MolarWeight=0.044,
            DfH=-393500,
            DfG_25degC=-394400,
            References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "CO2(g)";

      package CarbonDioxide_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.044,
            DfH=-412900,
            DfG_25degC=-386200,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "CO2(aq)";

       package Carbonate_aqueous =
            Interfaces.BaseSubstance (
            MolarWeight=0.06001,
          z=-2,
          DfH=-676300,
            DfG_25degC=-676300-298.15*(-497.065),
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "CO3--(aq)";

        package Electrone_solid =
          Interfaces.BaseSubstance (
            MolarWeight=5.4857990946e-7,
          z=-1,
          DfH=0,
            DfG_25degC=0,
            References={"http://physics.nist.gov/cgi-bin/cuu/Value?mme","To solve standard electo-chemical cell potentials"}) "e-(s)";

         package Iron2_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.05585,
          z=2,
          DfH=-87860,
            DfG_25degC=-87860-298.15*(-9.93),
            References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
        "}) "Fe++(aq)";

       package Iron3_aqueous =
            Interfaces.BaseSubstance (
            MolarWeight=0.05585,
          z=3,
          DfH=-47700,
            DfG_25degC=-47700-298.15*(-124.77),
            References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
"}) "Fe+++(aq)";

     package Glucose_solid =
          Interfaces.BaseSubstance (
            MolarWeight=0.1806,
            DfH=-1274500,
            DfG_25degC=-1274500-298.15*(-1220.66),
            References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
"}) "Glu(s)";

      package Hydrogen_gas =
          Interfaces.BaseSubstance (
            MolarWeight=0.00201588,
            z=0,
            DfH=0,
            DfG_25degC=0,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"}) "H2(g)";

      package CarbonicAcid_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.062027,
            DfH=-699700,
            DfG_25degC=-699700-298.15*(-256.582),
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H2CO3(aq)";

      package Water_gas =
          Interfaces.BaseSubstance (
            MolarWeight=0.018015,
            DfH=-241830,
            DfG_25degC=-228590,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H2O(g)";

      package Water_liquid =
          Interfaces.BaseSubstance (
            MolarWeight=0.018015,
            DfH=-285830,
            DfG_25degC=-237190,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H2O(l)";

      package DihydrogenPhosphate_aqueous =
            Interfaces.BaseSubstance (
            MolarWeight=0.095,
          z=-1,
          DfH=-1302480,
            DfG_25degC=-1302480-298.15*(-561.395),
            References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "H2PO4-(aq)";

      package Hydronium_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.019022,
            z=1,
            DfH=-285840,
            DfG_25degC= -285840 - 298.15* (-163.17),
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H3O+(aq)";

       package PhosphoricAcid_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.095,
            DfH=-1288000,
            DfG_25degC=-1288000-298.15*(-496.4),
            References={"https://en.wikipedia.org/wiki/Phosphoric_acid",
            "https://www.researchgate.net/publication/6600409_Standard_thermodynamic_properties_of_H3PO4%28aq%29_over_a_wide_range_of_temperatures_and_pressures"})
        "H3PO4(aq)";

      package Proton_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.001007,
            z=1,
            DfH=0,
            DfG_25degC=0,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "H+(aq)";
                 // as hypothetical HA <-> H+ + A- simplification of H2O + HA <-> H3O+ + A-";

       package Bicarbonate_aqueous =
            Interfaces.BaseSubstance (
            MolarWeight=0.06102,
          z=-1,
          DfH=-691100,
            DfG_25degC=-691100-298.15*(-348.82),
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "HCO3-(aq)";

          package HydrogenPhosphate_aqueous =
            Interfaces.BaseSubstance (
            MolarWeight=0.095,
          z=-2,
          DfH=-1298700,
            DfG_25degC=-1298700-298.15*(-686.232),
            References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "HPO4--(aq)";

          package Potassium_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.0391,
          z=1,
          DfH=-251200,
            DfG_25degC=-251200-298.15*(103.97),
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "K+(aq)";

          package Magnesium_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.0243,
          z=2,
          DfH=-461960,
            DfG_25degC=-461960-298.15*(-19.99),
            References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
","http://www.vias.org/genchem/standard_enthalpies_table.html"}) "Mg++(aq)";

        package Sodium_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.02299,
          z=1,
          DfH=-239660,
            DfG_25degC=-239660-298.15*(74.49),
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "Na+(aq)";

      package Amonium_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.01804,
          z=1,
          DfH=-132800,
            DfG_25degC=-132800-298.15*(-178.77),
            References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf
"}) "NH4+(aq)";

      package Oxygen_gas =
          Interfaces.BaseSubstance (
            MolarWeight=0.032,
            DfH=0,
            DfG_25degC=0,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"}) "O2(g)";

      package Oxygen_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.032,
            DfH=-11700,
            DfG_25degC=16320,
            References={
            "http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.pdf",
            "https://books.google.cz/books?id=dr-VBAAAQBAJ&pg=PA156&lpg=PA156&dq=Gibbs+energy+of+formation++%22O2(aq)%22&source=bl&ots=09N5CxY7OD&sig=hbsTXQvX59vXBqHUjFVVIZQpHCA&hl=cs&sa=X&ei=sDQtVaeUMMaRsAHpzYHgAg&redir_esc=y#v=onepage&q=Gibbs%20energy%20of%20formation%20%20%22O2(aq)%22&f=false"})
        "O2(aq)";
      package Hydroxide_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.017006,
            z=-1,
            DfH=-229940,
            DfG_25degC=-157300,
            References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
        "OH-(aq)";

           package Phosphate_aqueous =
            Interfaces.BaseSubstance (
            MolarWeight=0.095,
          z=-3,
          DfH=-1284070,
            DfG_25degC=-1284070-298.15*(-866.946),
            References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "PO4---(aq)";

        package Sulphates_aqueous =
          Interfaces.BaseSubstance (
            MolarWeight=0.09607,
          z=-2,
          DfH=-907500,
            DfG_25degC=-907500-298.15*(-555.123),
            References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
        "SO4--(aq)";

        //Some organic molecules: https://www.e-education.psu.edu/drupal6/files/be497b/pdf/Bioenergetics_AppA.pdf
    end Substances;

    model SimpleReaction
       extends Modelica.Icons.Example;

      constant Real K = 2 "Dissociation constant of the reaction";

      constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Components.Substance A(
        amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-64,-8},{-44,12}})));

      Components.Reaction reaction annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
      Components.Substance B(
        substance( DfG_25degC=-R*T_25degC*log(K)),
        amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{42,-8},{62,12}})));
      Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    equation

      connect(A.port_a, reaction.substrates[1]) annotation (Line(
          points={{-54,2},{-10,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(reaction.products[1], B.port_a) annotation (Line(
          points={{10,2},{52,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(A.solution, solution.solution) annotation (Line(
          points={{-60,-8},{-60,-92},{0,-92},{0,-100}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(B.solution, solution.solution) annotation (Line(points={{46,-8},{
              46,-92},{0,-92},{0,-100}}, smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=2e-005),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}),
                        graphics),
        __Dymola_experimentSetupOutput);
    end SimpleReaction;

    model SimpleReaction2
       extends Modelica.Icons.Example;

      constant Real Kb(unit="kg/mol") = 2
        "Molality based dissociation constant of the reaction with one more reactant";

      constant Real Kx(unit="1") = Kb*55.508
        "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

      constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Components.Substance A(amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
      Components.Reaction reaction(nS=2)
        annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Components.Substance B(amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Components.Substance C(substance(DfG_25degC=-R*T_25degC*log(Kx)), amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{48,-8},{68,12}})));
      Components.Solution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    equation

      connect(reaction.products[1], C.port_a) annotation (Line(
          points={{24,2},{24,2},{58,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(A.solution, solution.solution) annotation (Line(
          points={{-30,2},{-30,-90},{0,-90},{0,-100}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(C.solution, solution.solution) annotation (Line(points={{52,-8},{66,-8},
              {66,-90},{0,-90},{0,-100}},        smooth=Smooth.None));
      connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},{-30,
              -90},{0,-90},{0,-100}},    smooth=Smooth.None));
      connect(A.solution, B.solution)
        annotation (Line(points={{-30,2},{-30,-24}}, smooth=Smooth.None));
      connect(B.port_a, reaction.substrates[1]) annotation (Line(
          points={{-24,-14},{-10,-14},{-10,1.5},{4,1.5}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(A.port_a, reaction.substrates[2]) annotation (Line(
          points={{-24,12},{-10,12},{-10,2.5},{4,2.5}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=5e-006),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics),
        __Dymola_experimentSetupOutput);
    end SimpleReaction2;

    model ExothermicReaction

       extends Modelica.Icons.Example;

      constant Modelica.SIunits.MolarEnergy ReactionEnthalpy=-500;

      Components.Substance A( amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-56,-8},{-36,12}})));
      Components.Reaction reaction
        annotation (Placement(transformation(extent={{-8,-8},{12,12}})));
      Components.Substance B( amountOfSubstance_start=0.1, substance(DfH=ReactionEnthalpy))
        annotation (Placement(transformation(extent={{44,-8},{64,12}})));
      Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    equation

      connect(A.port_a, reaction.substrates[1]) annotation (Line(
          points={{-46,2},{-8,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(reaction.products[1], B.port_a) annotation (Line(
          points={{12,2},{54,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(B.solution, solution.solution) annotation (Line(
          points={{48,-8},{48,-84},{0,-84},{0,-100}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(A.solution, solution.solution)
        annotation (Line(points={{-52,-8},{-52,-84},{0,-84},{0,-100}}, smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=2e-005),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                        graphics),
        __Dymola_experimentSetupOutput);
    end ExothermicReaction;

    model Henry
       extends Modelica.Icons.Example;
      Components.GasSolubility CO2_new
        annotation (Placement(transformation(extent={{-86,-20},{-66,0}})));
      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,

      Components.Substance CO2_l_plasma(redeclare package substance =
            Substances.CarbonDioxide_aqueous)
        annotation (Placement(transformation(extent={{-86,-60},{-66,-40}})));
      Components.GasSolubility O2_new
        annotation (Placement(transformation(extent={{8,-20},{28,0}})));

      Sources.AirSubstance O2_g_n1(
        redeclare package substance=Substances.Oxygen_gas,
        PartialPressure=12665.626804425,
        TotalPressure=101325.0144354)
        annotation (Placement(transformation(extent={{32,34},{52,54}})));
      Components.Substance O2_l_plasma(redeclare package substance =
            Substances.Oxygen_aqueous)
        annotation (Placement(transformation(extent={{10,-60},{30,-40}})));
      Components.GasSolubility CO2_new1
        annotation (Placement(transformation(extent={{-50,-20},{-30,0}})));

      Sources.AirSubstance CO2_g_n2(
        redeclare package substance=Substances.CarbonDioxide_gas,
        PartialPressure=5332.8954966,
        TotalPressure=101325.0144354)
        annotation (Placement(transformation(extent={{-80,34},{-60,54}})));

      Components.Substance CO2_l_erythrocyte(redeclare package substance =
            Substances.CarbonDioxide_aqueous)
        annotation (Placement(transformation(extent={{-50,-62},{-30,-42}})));
      Components.GasSolubility O2_new1
        annotation (Placement(transformation(extent={{40,-18},{60,2}})));

      Components.Substance O2_l_erythrocyte(redeclare package substance =
            Substances.Oxygen_aqueous)
        annotation (Placement(transformation(extent={{40,-60},{60,-40}})));
      Components.Solution plasma(amountOfSolution_start=52.3)
        annotation (Placement(transformation(extent={{-104,-100},{96,100}})));
      Components.Solution redCells(amountOfSolution_start=39.7)
        annotation (Placement(transformation(extent={{-92,-90},{108,110}})));
      Components.GasSolubility O2_new2(useWaterCorrection=false)
        annotation (Placement(transformation(extent={{74,-20},{94,0}})));
      Components.Substance O2_l_erythrocyte1(substance(DfG_25degC=-Modelica.Constants.R*298.15*log(0.0013*0.018)))
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
          points={{-76,-20},{-76,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(CO2_new1.liquid_port, CO2_l_erythrocyte.port_a) annotation (Line(
          points={{-40,-20},{-40,-52}},
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
          points={{18,-20},{18,-50},{20,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(O2_new1.liquid_port, O2_l_erythrocyte.port_a) annotation (Line(
          points={{50,-18},{50,-50}},
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
          points={{84,-20},{84,-50}},
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
        experiment(StopTime=1));
    end Henry;

    model MichaelisMenten "Basic enzyme kinetics"
      extends Modelica.Icons.Example;

      //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
      Sources.AmbientMoleFraction
                              P(             substance(DfG_25degC=-Modelica.Constants.R*298.15*50),
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

          Components.Substance ES(substance(DfG_25degC=-Modelica.Constants.R*298.15*log(2/Km)),
          amountOfSubstance_start=tE/2)
            annotation (Placement(transformation(extent={{-12,-10},{8,10}})));
          Components.Substance E(amountOfSubstance_start=tE/2)
            annotation (Placement(transformation(extent={{-10,38},{10,58}})));
      Components.Reaction chemicalReaction(nS=2, Tau=2*AmountOfSolution*Modelica.Constants.R
            *298.15*log(2)/Vmax)
        annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));

      Components.Reaction chemicalReaction1(nP=2, Tau=2*AmountOfSolution*Modelica.Constants.R
            *298.15*(50 - log(2))/Vmax)
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
          points={{-22,0},{-2,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(ES.port_a, chemicalReaction1.substrates[1]) annotation (Line(
          points={{-2,0},{24,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(E.port_a, chemicalReaction.substrates[2]) annotation (Line(
          points={{0,48},{-52,48},{-52,0.5},{-42,0.5}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(E.port_a, chemicalReaction1.products[2]) annotation (Line(
          points={{0,48},{54,48},{54,0.5},{44,0.5}},
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
</html>"),
        experiment,
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                        graphics),
        __Dymola_experimentSetupOutput);
    end MichaelisMenten;

    model StandardElectrochemicalCellPotential
     extends Modelica.Icons.Example;

      Sources.PureSubstance Ag(redeclare package substance =
            Chemical.Examples.Substances.Silver_solid)
        annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
      Sources.PureSubstance Cl(
                              redeclare package substance =
            Chemical.Examples.Substances.Chloride_aqueous)
                                       annotation (Placement(transformation(extent={{-8,-36},{-28,-16}})));
      Sources.PureSubstance AgCl(redeclare package substance =
            Chemical.Examples.Substances.SilverChloride_solid)
        annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
      Sources.AirSubstance H2(
        redeclare package substance = Chemical.Examples.Substances.Hydrogen_gas,
        PartialPressure=100000,
        TotalPressure=100000) annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Sources.PureSubstance H(
        redeclare package substance =
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
                           electrone(redeclare package substance =
            Chemical.Examples.Substances.Electrone_solid)
                                     annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
      Sources.PureElectricParticle
                           electrone1(redeclare package substance =
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
            graphics));
    end StandardElectrochemicalCellPotential;

    model ElectrochemicalCell
     extends Modelica.Icons.Example;

      Sources.PureSubstance Ag(redeclare package substance =
            Chemical.Examples.Substances.Silver_solid)
        annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
      Components.Substance Cl(redeclare package substance =
            Chemical.Examples.Substances.Chloride_aqueous,
          amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-28,-36},{-8,-16}})));
      Sources.PureSubstance AgCl(redeclare package substance =
            Chemical.Examples.Substances.SilverChloride_solid)
        annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
      Sources.AirSubstance H2(
        redeclare package substance = Chemical.Examples.Substances.Hydrogen_gas,
        PartialPressure=100000,
        TotalPressure=100000) annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Components.Substance H(
        redeclare package substance =
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
                           electrone(redeclare package substance =
            Chemical.Examples.Substances.Electrone_solid)
                                     annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
      Sources.PureElectricParticle
                           electrone1(redeclare package substance =
            Chemical.Examples.Substances.Electrone_solid)
                                      annotation (Placement(transformation(extent={{86,-26},{66,-6}})));
    equation
      connect(Ag.port_a, electrodeReaction1.substrates[1]) annotation (Line(
          points={{-60,-18},{-42,-18},{-42,-4},{-40.5,-4}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Cl.port_a, electrodeReaction1.substrates[2]) annotation (Line(
          points={{-18,-26},{-39.5,-26},{-39.5,-4}},
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
          points={{20,-26},{52.5,-26},{52.5,-4}},
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
        constant Real KRx=KR*AmountOfSolutionIn1L "Mole fraction based KR";

      //Relative Gibbs formation energies of the substances in the system:
        constant Real RT = Modelica.Constants.R*T;
        constant Modelica.SIunits.MolarEnergy
          GO2aq=-RT*log(0.0013*0.018),
          GR0=0,                            GT0=GR0 -RT*log(L),
          GR1=GR0+GO2aq +RT*log(KRx/4),     GT1=GR1 -RT*log(c*L),
          GR2=GR1+GO2aq +RT*log(KRx/(3/2)), GT2=GR2 -RT*log(c^2*L),
          GR3=GR2+GO2aq +RT*log(KRx/(2/3)), GT3=GR3 -RT*log(c^3*L),
          GR4=GR3+GO2aq +RT*log(KRx*4),     GT4=GR4 -RT*log(c^4*L);

        Components.Substance oxygen_unbound(substance( DfG_25degC=GO2aq),
            amountOfSubstance_start(displayUnit="mol") = 1e-5)
          annotation (Placement(transformation(extent={{-56,-44},{-36,-24}})));

        Components.Substance T0(substance( DfG_25degC=GT0), amountOfSubstance_start=
              THb)
          annotation (Placement(transformation(extent={{34,78},{54,98}})));

        Components.Substance T1(substance( DfG_25degC=GT1),
            amountOfSubstance_start=1e-10)
          annotation (Placement(transformation(extent={{34,36},{54,56}})));

        Components.Substance T2(substance( DfG_25degC=GT2),
            amountOfSubstance_start=1e-18)
          annotation (Placement(transformation(extent={{34,-10},{54,10}})));

        Components.Substance R1(substance( DfG_25degC=GR1),
            amountOfSubstance_start=1e-15)
          annotation (Placement(transformation(extent={{-20,36},{0,56}})));

        Components.Substance R2(substance( DfG_25degC=GR2),
            amountOfSubstance_start=1e-20)
          annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

        Components.Substance T3(substance( DfG_25degC=GT3),
            amountOfSubstance_start=1e-25)
          annotation (Placement(transformation(extent={{34,-54},{54,-34}})));

        Components.Substance R3(substance( DfG_25degC=GR3),
            amountOfSubstance_start=1e-25)
          annotation (Placement(transformation(extent={{-20,-54},{0,-34}})));

        Components.Substance T4(substance( DfG_25degC=GT4),
            amountOfSubstance_start=1e-33)
          annotation (Placement(transformation(extent={{34,-92},{54,-72}})));

        Components.Substance R4(substance( DfG_25degC=GR4),
            amountOfSubstance_start=1e-31)
          annotation (Placement(transformation(extent={{-20,-92},{0,-72}})));

        Components.Substance R0(substance( DfG_25degC=GR0),
            amountOfSubstance_start=1e-10)
          annotation (Placement(transformation(extent={{-20,78},{0,98}})));

        Components.Reaction quaternaryForm annotation (Placement(transformation(extent={{4,78},{24,98}})));
        Components.Reaction oxyR1(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,64})));
        Components.Reaction oxyT1(nP=2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,64})));
        Components.Reaction oxyR2(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,22})));
        Components.Reaction oxyR3(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-24})));
        Components.Reaction oxyR4(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-66})));
        Components.Reaction oxyT2(nP=2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,22})));
        Components.Reaction oxyT3(nP=2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-24})));
        Components.Reaction oxyT4(nP=2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-66})));
        Components.Reaction quaternaryForm1 annotation (Placement(transformation(extent={{8,36},{28,56}})));
        Components.Reaction quaternaryForm2 annotation (Placement(transformation(extent={{8,-10},{28,10}})));
        Components.Reaction quaternaryForm3 annotation (Placement(transformation(extent={{8,-54},{28,-34}})));
        Components.Reaction quaternaryForm4 annotation (Placement(transformation(extent={{10,-92},{30,-72}})));

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
          temperature=T,
          PartialPressure(displayUnit="kPa") = 12000,
          TotalPressure(displayUnit="kPa") = 101325.0144354,
          redeclare package substance = Chemical.Examples.Substances.Oxygen_gas,
          usePartialPressureInput=true) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,22})));
        Components.GasSolubility gasSolubility(useWaterCorrection=false)
          annotation (Placement(transformation(extent={{-94,-16},{-74,4}})));
        Components.Solution solution(amountOfSolution_start=
              AmountOfSolutionIn1L, Isothermal=false)
          annotation (Placement(transformation(extent={{-56,-102},{100,104}})));

      equation
       //  sO2 = (R1.amountOfSubstance + 2*R2.amountOfSubstance + 3*R3.amountOfSubstance + 4*R4.amountOfSubstance + T1.amountOfSubstance + 2*T2.amountOfSubstance + 3*T3.amountOfSubstance + 4*T4.amountOfSubstance)/(4*totalAmountOfHemoglobin);
      //   totalAmountOfRforms = R0.amountOfSubstance + R1.amountOfSubstance + R2.amountOfSubstance + R3.amountOfSubstance + R4.amountOfSubstance;
      //   totalAmountOfTforms = T0.amountOfSubstance + T1.amountOfSubstance + T2.amountOfSubstance + T3.amountOfSubstance + T4.amountOfSubstance;

      //   totalAmountOfHemoglobin*normalizedState[1] = totalAmountOfRforms + totalAmountOfTforms;

        connect(quaternaryForm.products[1],T0. port_a) annotation (Line(
            points={{24,88},{44,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR1.substrates[1],R1. port_a) annotation (Line(
            points={{-10,54},{-10,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R1.port_a,oxyR2. products[1]) annotation (Line(
            points={{-10,46},{-10,32},{-10.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR2.substrates[1],R2. port_a) annotation (Line(
            points={{-10,12},{-10,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.substrates[1],R3. port_a) annotation (Line(
            points={{-10,-34},{-10,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.products[1],R2. port_a) annotation (Line(
            points={{-10.5,-14},{-10.5,-7},{-10,-7},{-10,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R3.port_a,oxyR4. products[1]) annotation (Line(
            points={{-10,-44},{-10,-56},{-10.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR4.substrates[1],R4. port_a) annotation (Line(
            points={{-10,-76},{-10,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT1.products[1],T0. port_a) annotation (Line(
            points={{44.5,74},{44.5,88},{44,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT1.substrates[1],T1. port_a) annotation (Line(
            points={{44,54},{44,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(T1.port_a,oxyT2. products[1]) annotation (Line(
            points={{44,46},{44,32},{44.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT3.substrates[1],T3. port_a) annotation (Line(
            points={{44,-34},{44,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(T3.port_a,oxyT4. products[1]) annotation (Line(
            points={{44,-44},{44,-56},{44.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT4.substrates[1],T4. port_a) annotation (Line(
            points={{44,-76},{44,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,quaternaryForm. substrates[1]) annotation (Line(
            points={{-10,88},{4,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,oxyR1. products[1]) annotation (Line(
            points={{-10,88},{-10,74},{-10.5,74}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R1.port_a,quaternaryForm1. substrates[1]) annotation (Line(
            points={{-10,46},{8,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm1.products[1],T1. port_a) annotation (Line(
            points={{28,46},{44,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R2.port_a,quaternaryForm2. substrates[1]) annotation (Line(
            points={{-10,0},{8,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R3.port_a,quaternaryForm3. substrates[1]) annotation (Line(
            points={{-10,-44},{8,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm3.products[1],T3. port_a) annotation (Line(
            points={{28,-44},{44,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R4.port_a,quaternaryForm4. substrates[1]) annotation (Line(
            points={{-10,-82},{10,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm4.products[1],T4. port_a) annotation (Line(
            points={{30,-82},{44,-82}},
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
            points={{-9.5,74},{-46,74},{-46,-34}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR2.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,32},{-46,32},{-46,-34}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,-14},{-46,-14},{-46,-34}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR4.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,-56},{-46,-56},{-46,-34}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT1.products[2])
                                            annotation (Line(
            points={{-46,-34},{-46,74},{43.5,74}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT2.products[2])
                                            annotation (Line(
            points={{-46,-34},{-46,32},{43.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT3.products[2])
                                            annotation (Line(
            points={{-46,-34},{-46,-14},{43.5,-14}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT4.products[2])
                                            annotation (Line(
            points={{-46,-34},{-46,-56},{43.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(O2_in_air.port_a, gasSolubility.gas_port) annotation (Line(
            points={{-84,12},{-84,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(gasSolubility.liquid_port, oxygen_unbound.port_a) annotation (Line(
            points={{-84,-16},{-84,-34},{-46,-34}},
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
            points={{28,0},{44,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT2.substrates[1], T2.port_a) annotation (Line(
            points={{44,12},{44,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(T2.port_a, oxyT3.products[1]) annotation (Line(
            points={{44,0},{44,-14},{44.5,-14}},
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
            StopTime=1e-005,
            Tolerance=1e-005,
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
<p><i>2013</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics),
          __Dymola_experimentSetupOutput);
      end Allosteric_Hemoglobin_MWC;

      model Allosteric_Hemoglobin2_MWC
        "Monod,Wyman,Changeux (1965) - The same allosteric hemoglobin model as Allosteric_Hemoglobin_MWC implemented by Speciation blocks"

       extends Modelica.Icons.Example;

        parameter Physiolibrary.Types.MolarEnergy DfHT=10000
          "Enthalpy of formation of heme oxygenation in T hemoglobin form";
        parameter Physiolibrary.Types.MolarEnergy DfHR=20000
          "Enthalpy of formation of heme oxygenation in R hemoglobin form";
        parameter Physiolibrary.Types.MolarEnergy DfHL=-1000
          "Enthalpy of formation of reaction T->R as hemoglobin tetramer structure change";

        parameter Physiolibrary.Types.Fraction L=7.0529*10^6
          "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
        parameter Physiolibrary.Types.Fraction c=0.00431555
          "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
        parameter Physiolibrary.Types.Concentration KR=0.000671946
          "oxygen dissociation on relaxed(R) hemoglobin subunit";
                                                                    //*7.875647668393782383419689119171e-5
                                                                  //10.500001495896 7.8756465463794e-05

        parameter Physiolibrary.Types.Concentration KT=KR/c
          "oxygen dissociation on tensed(T) hemoglobin subunit";

        parameter Physiolibrary.Types.AmountOfSubstance totalAmountOfHemoglobin=1;

        Components.Reaction quaternaryForm(
          K=L,
          TK=310.15,
          DfH=DfHL)
          annotation (Placement(transformation(extent={{-2,-76},{18,-56}})));
        Components.Speciation R0_in_R(NumberOfSubunits=4, useInternalHeatsInput=
             true)
          annotation (Placement(transformation(extent={{-30,-68},{-50,-48}})));
        Components.Speciation T0_in_T(NumberOfSubunits=4, useInternalHeatsInput=
             true)
          annotation (Placement(transformation(extent={{70,-66},{50,-46}})));
        Components.Substance OxyRHm[4](
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          each isDependent=true,
          each solute_start=4e-19,
          each DfH=-DfHL/4 - DfHR)
          "Oxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-96,-18},{-76,2}})));
        Components.Reaction oxygenation_R[4](
          each K=KR,
          each nP=2,
          each TK=310.15,
          each DfH=DfHR)
          annotation (Placement(transformation(extent={{-68,-18},{-48,2}})));
        Components.Substance DeoxyRHm[4](each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          each solute_start=4e-11,
          each DfH=-DfHL/4)
          "Deoxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-40,-18},{-20,2}})));
        Components.Substance OxyTHm[4](
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isDependent={false,true,true,true},
          each DfH=-DfHT,
          each solute_start=1e-14)
          "Oxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{14,-18},{34,2}})));
        Components.Reaction oxygenation_T[4](
          each K=KT,
          each nP=2,
          each DfH=DfHT,
          each TK=310.15)
          annotation (Placement(transformation(extent={{42,-18},{62,2}})));
        Components.Substance DeoxyTHm[4](        each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          each solute_start=0.00025,
          each DfH=0)
          "Deoxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{70,-18},{90,2}})));

        Components.Substance oxygen_unbound(
                                           Simulation=Physiolibrary.Types.SimulationType.SteadyState,
                                                                                        solute_start=0.000001
              *7.875647668393782383419689119171e-5)
          annotation (Placement(transformation(extent={{-2,6},{18,26}})));
        Modelica.Blocks.Sources.Clock clock(offset=10)
          annotation (Placement(transformation(extent={{-40,74},{-20,94}})));
        Modelica.Blocks.Math.Add add[4] annotation (Placement(transformation(
              extent={{-4,-4},{4,4}},
              rotation=270,
              origin={-58,-36})));
        Modelica.Blocks.Math.Add add1[4] annotation (Placement(transformation(
              extent={{-4,-4},{4,4}},
              rotation=270,
              origin={30,-48})));
        Physiolibrary.Chemical.Sources.UnlimitedGasStorage oxygen_in_air(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          usePartialPressureInput=true,
          isIsolatedInSteadyState=false,
          T=310.15) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={8,68})));
        Physiolibrary.Chemical.Components.GasSolubility partialPressure1(
          kH_T0(displayUnit="(mmol/l)/kPa at 25degC") = 0.026029047188736,
          T=310.15,
          C=1700) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                origin={8,40})));
        Physiolibrary.SteadyStates.Components.MolarConservationLaw totalHb(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          Total(displayUnit="mol") = totalAmountOfHemoglobin,
          n=2)
          annotation (Placement(transformation(extent={{72,-84},{92,-64}})));
        Modelica.Blocks.Math.Sum sum1(nin=8, k=(1/4)*ones(8))
                                             annotation (Placement(transformation(
              extent={{-4,-4},{4,4}},
              rotation=270,
              origin={-72,-74})));
        Modelica.Blocks.Math.Division sO2_ "hemoglobin oxygen saturation"
          annotation (Placement(transformation(extent={{-62,-88},{-52,-78}})));
        Modelica.Blocks.Math.Sum internalHeat(nin=2) "hemoglobin enthalpy heat"
          annotation (Placement(transformation(
              extent={{-4,-4},{4,4}},
              origin={8,-90})));
        Modelica.Blocks.Math.Add add2[
                                     4] annotation (Placement(transformation(
              extent={{-5,-5},{5,5}},
              rotation=270,
              origin={-73,-37})));
        Modelica.Blocks.Math.Add add3[4] annotation (Placement(transformation(
              extent={{-5,-5},{5,5}},
              rotation=270,
              origin={47,-39})));
      equation

        connect(R0_in_R.specificForm, quaternaryForm.substrates[1])
                                                         annotation (Line(
            points={{-50,-66},{-2,-66}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm.products[1], T0_in_T.specificForm)
                                                       annotation (Line(
            points={{18,-66},{34,-66},{34,-64},{50,-64}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(OxyTHm.q_out, oxygenation_T.substrates[1])
                                                 annotation (Line(
            points={{24,-8},{42,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_T.products[1], DeoxyTHm.q_out)
                                               annotation (Line(
            points={{62,-8.5},{72,-8.5},{72,-8},{80,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(OxyTHm.solute, add1.u2) annotation (Line(
            points={{30,-18},{30,-24},{27.6,-24},{27.6,-43.2}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(add1.u1, DeoxyTHm.solute) annotation (Line(
            points={{32.4,-43.2},{32.4,-24},{86,-24},{86,-18}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(partialPressure1.q_out, oxygen_in_air.q_out)
                                                  annotation (Line(
            points={{8,50},{8,58}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(partialPressure1.q_in, oxygen_unbound.q_out) annotation (Line(
            points={{8,32},{8,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(clock.y, oxygen_in_air.partialPressure) annotation (Line(
            points={{-19,84},{8,84},{8,78}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(add.y, R0_in_R.amountOfSubunit) annotation (Line(
            points={{-58,-40.4},{-58,-58},{-48,-58}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(OxyRHm.solute, add.u2) annotation (Line(
            points={{-80,-18},{-80,-24},{-60.4,-24},{-60.4,-31.2}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(DeoxyRHm.solute, add.u1) annotation (Line(
            points={{-24,-18},{-30,-18},{-30,-24},{-55.6,-24},{-55.6,-31.2}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(OxyRHm.q_out, oxygenation_R.substrates[1]) annotation (Line(
            points={{-86,-8},{-68,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(DeoxyRHm.q_out, R0_in_R.specificSubunitForm) annotation (Line(
            points={{-30,-8},{-40,-8},{-40,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_R.products[1], DeoxyRHm.q_out) annotation (Line(
            points={{-48,-8.5},{-40,-8.5},{-40,-8},{-30,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_R[1].products[2], oxygen_unbound.q_out) annotation (Line(
            points={{-48,-7.5},{-34,-7.5},{-34,16},{8,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_R[2].products[2], oxygen_unbound.q_out) annotation (Line(
            points={{-48,-7.5},{-34,-7.5},{-34,16},{8,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_R[3].products[2], oxygen_unbound.q_out) annotation (Line(
            points={{-48,-7.5},{-34,-7.5},{-34,16},{8,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_R[4].products[2], oxygen_unbound.q_out) annotation (Line(
            points={{-48,-7.5},{-34,-7.5},{-34,16},{8,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_T[1].products[2], oxygen_unbound.q_out) annotation (Line(
            points={{62,-7.5},{78,-7.5},{78,16},{8,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_T[2].products[2], oxygen_unbound.q_out) annotation (Line(
            points={{62,-7.5},{78,-7.5},{78,16},{8,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
       connect(oxygenation_T[3].products[2], oxygen_unbound.q_out) annotation (Line(
            points={{62,-7.5},{78,-7.5},{78,16},{8,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
            connect(oxygenation_T[4].products[2], oxygen_unbound.q_out) annotation (Line(
            points={{62,-7.5},{78,-7.5},{78,16},{8,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(T0_in_T.specificSubunitForm, DeoxyTHm.q_out)
                                                     annotation (Line(
            points={{60,-46},{84,-46},{84,-8},{80,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(add1.y, T0_in_T.amountOfSubunit) annotation (Line(
            points={{30,-52.4},{30,-56},{52,-56}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(R0_in_R.amount, totalHb.fragment[1]) annotation (Line(
            points={{-40,-66},{-40,-79},{72,-79}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(T0_in_T.amount, totalHb.fragment[2]) annotation (Line(
            points={{60,-64},{60,-64},{60,-77},{72,-77}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(OxyRHm.solute, sum1.u[1:4]) annotation (Line(
            points={{-80,-18},{-86,-18},{-86,-62},{-72,-62},{-72,-69.2},{-72.1,
                -69.2}},
            color={0,0,127},
            smooth=Smooth.Bezier));

        connect(OxyTHm.solute, sum1.u[5:8]) annotation (Line(
            points={{30,-18},{30,-60},{-71.3,-60},{-71.3,-69.2}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(sO2_.u1, sum1.y) annotation (Line(
            points={{-63,-80},{-72,-80},{-72,-78.4}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(totalHb.totalAmountOfSubstance, sO2_.u2) annotation (Line(
            points={{92,-78},{100,-78},{100,-100},{-76,-100},{-76,-86},{-63,-86}},
            color={0,0,127},
            smooth=Smooth.Bezier));
        connect(R0_in_R.internalHeat, internalHeat.u[1]) annotation (Line(
            points={{-34,-66},{-34,-90.4},{3.2,-90.4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(T0_in_T.internalHeat, internalHeat.u[2]) annotation (Line(
            points={{66,-64},{66,-74},{-24,-74},{-24,-89.6},{3.2,-89.6}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(add3.y, T0_in_T.subunitInternalHeat) annotation (Line(
            points={{47,-44.5},{47,-50},{52,-50}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(add2.y, R0_in_R.subunitInternalHeat) annotation (Line(
            points={{-73,-42.5},{-73,-52},{-48,-52}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(add3.u1, DeoxyTHm.internalHeat) annotation (Line(
            points={{50,-33},{50,-30},{91.6,-30},{91.6,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(OxyTHm.internalHeat, add3.u2) annotation (Line(
            points={{35.6,-12},{35.6,-30},{44,-30},{44,-33}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(OxyRHm.internalHeat, add2.u2) annotation (Line(
            points={{-74.4,-12},{-74.4,-28},{-76,-28},{-76,-31}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(DeoxyRHm.internalHeat, add2.u1) annotation (Line(
            points={{-18.4,-12},{-18.4,-28},{-70,-28},{-70,-31}},
            color={0,0,127},
            smooth=Smooth.None));
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
          experiment(
            StopTime=15000,
            Tolerance=1e-014,
            __Dymola_Algorithm="Euler"),
          Documentation(revisions=
                        "<html>
<p><i>2013</i></p>
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
</html>"));
      end Allosteric_Hemoglobin2_MWC;

      model Hemoglobin_MKM_Specie "Part of model Hemoglobin_MKM_Adair"

      parameter Boolean loadStarts
          "Start values of state variables from data file (to help with initialization)";
      parameter Boolean storeState
          "Save state variables at the end of simulation";
      constant String dirName = Modelica.Utilities.Files.loadResource("modelica://Physiolibrary/Resources/Data/Hemoglobin_MKM")
          "Directory to load start gues values and store final simulation values";

      parameter Real[4] pKz
          "Dissociation coefficient of reaction z (Val1 amino terminal protonation)";
      parameter Real[4] pKc
          "Dissociation coefficient of reaction c (Val1 amino terminal carbamination)";
      parameter Real[4] pKh
          "Dissociation coefficient of reaction h (other Bohr protonation reactions of side chains)";

      parameter Physiolibrary.Types.MolarEnergy[4] dH_HbuANH2
          "Standard enthalpy of deprotonated and decarboxylated hemoglobin subunit";
      parameter Physiolibrary.Types.MolarEnergy[4] dHz
          "Enthalpy of reaction z (Val1 amino terminal protonation)";
      parameter Physiolibrary.Types.MolarEnergy[4] dHc
          "Enthalpy of reaction c (Val1 amino terminal carbamination)";
      parameter Physiolibrary.Types.MolarEnergy[4] dHh
          "Enthalpy of reaction h (other Bohr protonation reactions of side chains)";

      parameter Boolean isDependent=false
          "contains dependent equation (if solver is not smart enough)";

      Physiolibrary.Chemical.Interfaces.ChemicalPort_a Hbtn
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
          Physiolibrary.Chemical.Components.Substance Hbu_A_NH3[4](
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          dH=dH_HbuANH2 - dHz,
          each dirName=dirName,
          each LOAD_STARTS=loadStarts,
          each SAVE_RESULTS=storeState,
          each solute_start=1e-06)
          annotation (Placement(transformation(extent={{-32,70},{-12,90}})));
      Physiolibrary.Chemical.Components.Substance Hbu_AH_NH3[4](
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          each dirName=dirName,
          each LOAD_STARTS=loadStarts,
          each SAVE_RESULTS=storeState,
          each solute_start=1e-06,
          dH=dH_HbuANH2 - dHh - dHz)
          annotation (Placement(transformation(extent={{54,70},{74,90}})));
      Physiolibrary.Chemical.Components.Substance Hbu_A_NH2[4](
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isDependent={isDependent,true,true,true},
          each dirName=dirName,
          each LOAD_STARTS=loadStarts,
          each SAVE_RESULTS=storeState,
          each solute_start=1e-06,
          dH=dH_HbuANH2)
          annotation (Placement(transformation(extent={{-32,-2},{-12,18}})));
      Physiolibrary.Chemical.Components.Substance Hbu_AH_NH2[4](
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          each dirName=dirName,
          each LOAD_STARTS=loadStarts,
          each SAVE_RESULTS=storeState,
          each solute_start=1e-06,
          dH=dH_HbuANH2 - dHh)
          annotation (Placement(transformation(extent={{54,-2},{74,18}})));
      Physiolibrary.Chemical.Components.Substance Hbu_A_NHCOO[4](
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          dH=dH_HbuANH2 + dHc,
          each dirName=dirName,
          each LOAD_STARTS=loadStarts,
          each SAVE_RESULTS=storeState,
          each solute_start=1e-06)
          annotation (Placement(transformation(extent={{-32,-84},{-12,-64}})));
      Physiolibrary.Chemical.Components.Substance Hbu_AH_NHCOO[4](
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          each dirName=dirName,
          each LOAD_STARTS=loadStarts,
          each SAVE_RESULTS=storeState,
          dH=dH_HbuANH2 + dHc,
          each solute_start=1e-06)
          annotation (Placement(transformation(extent={{54,-84},{74,-64}})));
      Physiolibrary.Chemical.Components.ChemicalReaction h2[4](
          each nS=1,
          each nP=2,
          K=fill(10, 4) .^ (-pKh .+ 3),
          each TK=310.15,
          dH=dHh)
          annotation (Placement(transformation(extent={{32,-2},{12,18}})));
      Physiolibrary.Chemical.Components.ChemicalReaction z1[4](
          each nP=2,
          K=fill(10, 4) .^ (-pKz .+ 3),
          dH=dHz,
          each TK=310.15) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-22,44})));
      Physiolibrary.Chemical.Components.ChemicalReaction z2[4](
          each nP=2,
          K=fill(10, 4) .^ (-pKz .+ 3),
          each TK=310.15,
          dH=dHz) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={64,44})));
      Physiolibrary.Chemical.Components.ChemicalReaction c1[4](
          each nS=2,
          each nP=2,
          K=fill(10, 4) .^ (-pKc .+ 3),
          each TK=310.15,
          dH=dHc) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-22,-34})));
      Physiolibrary.Chemical.Components.ChemicalReaction c2[4](
          each nS=2,
          each nP=2,
          K=fill(10, 4) .^ (-pKc .+ 3),
          each TK=310.15,
          dH=dHc) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={64,-34})));
      Modelica.Blocks.Math.Sum totalAmounts[4](each nin=6) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-64,62})));
        Physiolibrary.Chemical.Interfaces.ChemicalPort_a H(conc(nominal=10^(-7.2
                 + 3))) "hydrogen ions"
          annotation (Placement(transformation(extent={{90,76},{110,96}})));
        Physiolibrary.Chemical.Interfaces.ChemicalPort_a CO2
          annotation (Placement(transformation(extent={{90,-70},{110,-50}})));
        Physiolibrary.Chemical.Components.Speciation Hb_tn(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          NumberOfSubunits=4,
          useInternalHeatsInput=true)
          annotation (Placement(transformation(extent={{-54,-22},{-74,-2}})));
      Physiolibrary.Types.RealIO.AmountOfSubstanceOutput tHb_u annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-100,-78})));
      Physiolibrary.Types.RealIO.EnergyOutput internalHeat "internal heat"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-100,-98})));
      Modelica.Blocks.Math.Sum totalHeats[4](each nin=6) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-64,32})));
      equation
      connect(Hbu_AH_NH3.q_out, z2.substrates[1]) annotation (Line(
          points={{64,80},{64,54}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(Hbu_A_NH3.q_out, z1.substrates[1]) annotation (Line(
          points={{-22,80},{-22,54}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(z1.products[1], Hbu_A_NH2.q_out) annotation (Line(
          points={{-22.5,34},{-22.5,22},{-22,22},{-22,8}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(z2.products[1], Hbu_AH_NH2.q_out) annotation (Line(
          points={{63.5,34},{63.5,22},{64,22},{64,8}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(h2.substrates[1], Hbu_AH_NH2.q_out) annotation (Line(
          points={{32,8},{64,8}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(Hbu_A_NH2.q_out, c1.substrates[1]) annotation (Line(
          points={{-22,8},{-22,-24},{-22.5,-24}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(Hbu_AH_NH2.q_out, c2.substrates[1]) annotation (Line(
          points={{64,8},{64,-24},{63.5,-24}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(c1.products[1], Hbu_A_NHCOO.q_out) annotation (Line(
          points={{-22.5,-44},{-22.5,-60},{-22,-60},{-22,-74}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(c2.products[1], Hbu_AH_NHCOO.q_out) annotation (Line(
          points={{63.5,-44},{63.5,-60},{64,-60},{64,-74}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
        connect(Hbu_A_NH3.solute, totalAmounts.u[1]) annotation (Line(
            points={{-16,70},{-44,70},{-44,63.6667},{-52,63.6667}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_AH_NH3.solute, totalAmounts.u[2]) annotation (Line(
            points={{70,70},{-2,70},{-2,63},{-52,63}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_A_NH2.solute, totalAmounts.u[3]) annotation (Line(
            points={{-16,-2},{-44,-2},{-44,62.3333},{-52,62.3333}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_AH_NH2.solute, totalAmounts.u[4]) annotation (Line(
            points={{70,-2},{-2,-2},{-2,61.6667},{-52,61.6667}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_A_NHCOO.solute, totalAmounts.u[5]) annotation (Line(
            points={{-16,-84},{-44,-84},{-44,61},{-52,61}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_AH_NHCOO.solute, totalAmounts.u[6]) annotation (Line(
            points={{70,-84},{-2,-84},{-2,60.3333},{-52,60.3333}},
            color={0,0,127},
            smooth=Smooth.None));

      connect(Hbu_A_NH2.q_out, h2.products[1]) annotation (Line(
          points={{-22,8},{-10,8},{-10,7.5},{12,7.5}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));

        connect(Hb_tn.specificForm, Hbtn) annotation (Line(
            points={{-74,-20},{-86,-20},{-86,0},{-98,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(totalAmounts.y, Hb_tn.amountOfSubunit) annotation (Line(
            points={{-75,62},{-78,62},{-78,-12},{-72,-12}},
            color={0,0,127},
            smooth=Smooth.None));

        for i in 1:4 loop
          connect(z1[i].products[2], H) annotation (Line(
            points={{-21.5,34},{-21.5,28},{-4,28},{-4,96},{88,96},{88,86},{100,
                86}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(z2[i].products[2], H) annotation (Line(
            points={{64.5,34},{64.5,28},{-4,28},{-4,96},{88,96},{88,86},{100,86}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(c1[i].products[2], H) annotation (Line(
            points={{-21.5,-44},{-21.5,-50},{-4,-50},{-4,96},{88,96},{88,86},{
                100,86}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(c2[i].products[2], H) annotation (Line(
            points={{64.5,-44},{64.5,-50},{-4,-50},{-4,96},{88,96},{88,86},{100,
                86}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

          connect(H, h2[i].products[2]) annotation (Line(
            points={{100,86},{88,86},{88,96},{-4,96},{-4,8.5},{12,8.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

          connect(CO2, c2[i].substrates[2]) annotation (Line(
            points={{100,-60},{88,-60},{88,-20},{64.5,-20},{64.5,-24}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2, c1[i].substrates[2]) annotation (Line(
            points={{100,-60},{88,-60},{88,-20},{-21.5,-20},{-21.5,-24}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        end for;
        connect(Hb_tn.specificSubunitForm, Hbu_A_NH2.q_out) annotation (Line(
            points={{-64,-2},{-64,8},{-22,8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(tHb_u, Hb_tn.amount) annotation (Line(
            points={{-100,-78},{-64,-78},{-64,-20}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(Hb_tn.internalHeat, internalHeat) annotation (Line(
            points={{-58,-20},{-58,-98},{-100,-98}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(Hbu_A_NH3.internalHeat, totalHeats.u[1]) annotation (Line(
            points={{-10.4,76},{-44,76},{-44,33.6667},{-52,33.6667}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_AH_NH3.internalHeat, totalHeats.u[2]) annotation (Line(
            points={{75.6,76},{-2,76},{-2,33},{-52,33}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_A_NH2.internalHeat, totalHeats.u[3]) annotation (Line(
            points={{-10.4,4},{-44,4},{-44,32.3333},{-52,32.3333}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_AH_NH2.internalHeat, totalHeats.u[4]) annotation (Line(
            points={{75.6,4},{-2,4},{-2,31.6667},{-52,31.6667}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_A_NHCOO.internalHeat, totalHeats.u[5]) annotation (Line(
            points={{-10.4,-78},{-44,-78},{-44,31},{-52,31}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hbu_AH_NHCOO.internalHeat, totalHeats.u[6]) annotation (Line(
            points={{75.6,-78},{-2,-78},{-2,30.3333},{-52,30.3333}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(totalHeats.y, Hb_tn.subunitInternalHeat) annotation (Line(
            points={{-75,32},{-76,32},{-76,-6},{-72,-6}},
            color={0,0,127},
            smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>

<p>[1] Mateják M, Kulhánek T, Matouaek S. Adair-Based Hemoglobin Equilibrium with Oxygen, Carbon Dioxide and Hydrogen Ion Activity. Scandinavian Journal of Clinical &AMP; Laboratory Investigation; 2015</p>

<p>[2] Bauer C, Schr&ouml;der E. Carbamino compounds of haemoglobin in human adult and foetal blood. The Journal of physiology 1972;227:457-71.</p>

<p>[3] Siggaard-Andersen O. Oxygen-Linked Hydrogen Ion Binding of Human Hemoglobin. Effects of Carbon Dioxide and 2, 3-Diphosphoglycerate I. Studies on Erythrolysate. Scandinavian Journal of Clinical &AMP; Laboratory Investigation 1971;27:351-60.</p>

</html>"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics));
      end Hemoglobin_MKM_Specie;

      model Hemoglobin_MKM_Adair "Matejak,Kulhanek,Matousek (2014)"
        extends Modelica.Icons.Example;

        constant Real pKzD=7.73,pKcD=7.54,pKhD=7.52;
        constant Real pKzO=7.25,pKcO=8.35,pKhO=6.89;
        constant Physiolibrary.Types.MolarEnergy dHzD=-51400;
        constant Physiolibrary.Types.MolarEnergy dHzO=7700;
        constant Physiolibrary.Types.MolarEnergy dHcD=59100;
        constant Physiolibrary.Types.MolarEnergy dHcO=-41100;
        constant Physiolibrary.Types.MolarEnergy dHhD=49000;
        constant Physiolibrary.Types.MolarEnergy dHhO=-105000;
        constant Physiolibrary.Types.MolarEnergy dHo=50000;
        constant Physiolibrary.Types.MolarEnergy dH_HbuDANH2=0;
        // dHhD=0, dHhO=-104000, dHo=12700, dH_HbuDANH2=0;                           // dHhD=48600, dHhO=-104000, dHo=50000, dH_HbuDANH2=0;

        parameter Boolean storeResults=false;
        parameter Boolean loadStarts=true;

        Physiolibrary.Chemical.Components.ChemicalReaction K1(
          K=0.0121,
          nS=1,
          nP=2) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-44,68})));
        Physiolibrary.Chemical.Components.ChemicalReaction K2(
          K=0.0117,
          nS=1,
          nP=2) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-46,28})));
        Physiolibrary.Chemical.Components.ChemicalReaction K3(
          K=0.0871,
          nS=1,
          nP=2) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-48,-18})));
        Physiolibrary.Chemical.Components.ChemicalReaction K4(
          K=0.000386,
          nS=1,
          nP=2) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-50,-60})));
        Hemoglobin_MKM_Specie Hb0(
          pKz=fill(pKzD, 4),
          pKc=fill(pKcD, 4),
          pKh=fill(pKhD, 4),
          isDependent=true,
          dH_HbuANH2(displayUnit="kJ/mol") = fill(dH_HbuDANH2, 4),
          dHz(displayUnit="kJ/mol") = fill(dHzD, 4),
          dHc(displayUnit="kJ/mol") = fill(dHcD, 4),
          dHh(displayUnit="kJ/mol") = fill(dHhD, 4),
          storeState=storeResults,
          loadStarts=loadStarts)
          annotation (Placement(transformation(extent={{-24,78},{-4,98}})));
        Hemoglobin_MKM_Specie Hb1(
          pKz=cat(  1,
                    fill(pKzD, 3),
                    fill(pKzO, 1)),
          pKc=cat(  1,
                    fill(pKcD, 3),
                    fill(pKcO, 1)),
          pKh=cat(  1,
                    fill(pKhD, 3),
                    fill(pKhO, 1)),
          dH_HbuANH2(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dH_HbuDANH2, 3),
                  fill(dH_HbuDANH2 - dHo, 1)),
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
          storeState=storeResults,
          loadStarts=loadStarts)
          annotation (Placement(transformation(extent={{-24,40},{-4,60}})));
        Hemoglobin_MKM_Specie Hb2(
          pKz=cat(  1,
                    fill(pKzD, 2),
                    fill(pKzO, 2)),
          pKc=cat(  1,
                    fill(pKcD, 2),
                    fill(pKcO, 2)),
          pKh=cat(  1,
                    fill(pKhD, 2),
                    fill(pKhO, 2)),
          dH_HbuANH2(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dH_HbuDANH2, 2),
                  fill(dH_HbuDANH2 - dHo, 2)),
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
          storeState=storeResults,
          loadStarts=loadStarts)
          annotation (Placement(transformation(extent={{-24,0},{-4,20}})));
        Hemoglobin_MKM_Specie Hb3(
          pKz=cat(  1,
                    fill(pKzD, 1),
                    fill(pKzO, 3)),
          pKc=cat(  1,
                    fill(pKcD, 1),
                    fill(pKcO, 3)),
          pKh=cat(  1,
                    fill(pKhD, 1),
                    fill(pKhO, 3)),
          dH_HbuANH2(displayUnit="kJ/mol") = cat(
                  1,
                  fill(dH_HbuDANH2, 1),
                  fill(dH_HbuDANH2 - dHo, 3)),
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
          storeState=storeResults,
          loadStarts=loadStarts)
          annotation (Placement(transformation(extent={{-24,-44},{-4,-24}})));
        Hemoglobin_MKM_Specie Hb4(
          pKz=fill(pKzO, 4),
          pKc=fill(pKcO, 4),
          pKh=fill(pKhO, 4),
          dH_HbuANH2(displayUnit="kJ/mol") = fill(dH_HbuDANH2 - dHo, 4),
          dHz(displayUnit="kJ/mol") = fill(dHzO, 4),
          dHc(displayUnit="kJ/mol") = fill(dHcO, 4),
          dHh(displayUnit="kJ/mol") = fill(dHhO, 4),
          storeState=storeResults,
          loadStarts=loadStarts)
          annotation (Placement(transformation(extent={{-24,-88},{-4,-68}})));
        Physiolibrary.Chemical.Sources.UnlimitedGasStorage CO2(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isIsolatedInSteadyState=false,
          PartialPressure=0)
          annotation (Placement(transformation(extent={{96,72},{76,92}})));
        Physiolibrary.Chemical.Sources.UnlimitedSolutionStorage pH(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isIsolatedInSteadyState=false,
          Conc(displayUnit="mol/l") = 10^(-7 + 3)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={34,82})));
        Physiolibrary.SteadyStates.Components.MolarConservationLaw totalHemoglobin(
          n=5,
          Total(displayUnit="mol") = 1,
          Simulation=Physiolibrary.Types.SimulationType.SteadyState)
          annotation (Placement(transformation(extent={{44,6},{64,26}})));
        Modelica.Blocks.Math.Sum sO2(nin=4, k={4/4,3/4,2/4,1/4})
          annotation (Placement(transformation(extent={{62,-30},{82,-10}})));
        Physiolibrary.Chemical.Components.Substance oxygen_unbound(Simulation=
              Physiolibrary.Types.SimulationType.SteadyState, solute_start=
              1e-08)
          annotation (Placement(transformation(extent={{-94,-28},{-74,-8}})));
        Modelica.Blocks.Sources.Clock clock(offset=10)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,70})));
        Physiolibrary.Chemical.Sources.UnlimitedGasStorage oxygen_in_air(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isIsolatedInSteadyState=false,
          PartialPressure(displayUnit="Pa") = 10,
          usePartialPressureInput=true,
          T=310.15) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,34})));
        Physiolibrary.Chemical.Components.GasSolubility partialPressure1(
          kH_T0(displayUnit="(mmol/l)/kPa at 25degC") = 0.026029047188736,
          C=1700,
          T=310.15) annotation (Placement(transformation(extent={{-10,-10},{10,
                  10}}, origin={-84,6})));
        Modelica.Blocks.Math.Sum internalHeat(nin=5)
          annotation (Placement(transformation(extent={{44,-62},{64,-42}})));
        Modelica.Blocks.Math.Gain gain(k=4)
          annotation (Placement(transformation(extent={{38,-92},{46,-84}})));
        Modelica.Blocks.Continuous.Der der1
          annotation (Placement(transformation(extent={{52,-80},{60,-72}})));
        Modelica.Blocks.Continuous.Der der2
          annotation (Placement(transformation(extent={{52,-92},{60,-84}})));
        Modelica.Blocks.Math.Division derInternalHeat_per_derO2
          annotation (Placement(transformation(extent={{72,-92},{92,-72}})));
        Physiolibrary.Chemical.Components.GasSolubility partialPressure2(
          T=310.15,
          kH_T0(displayUnit="(mmol/l)/kPa at 25degC") = 0.60734443440384,
          C=2400) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                origin={68,62})));
        Physiolibrary.Chemical.Components.Substance CO2_unbound(Simulation=
              Physiolibrary.Types.SimulationType.SteadyState, solute_start=
              0.0012)
          annotation (Placement(transformation(extent={{58,30},{78,50}})));
      equation
        connect(oxygen_unbound.q_out, K2.products[1]) annotation (Line(
            points={{-84,-18},{-62,-18},{-62,42},{-46,42},{-46,38},{-46.5,38}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.q_out, K3.products[1]) annotation (Line(
            points={{-84,-18},{-62,-18},{-62,0},{-48.5,0},{-48.5,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(K1.products[1], oxygen_unbound.q_out) annotation (Line(
            points={{-44.5,78},{-44.5,80},{-62,80},{-62,-18},{-84,-18}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.q_out, K4.products[1]) annotation (Line(
            points={{-84,-18},{-62,-18},{-62,-44},{-50.5,-44},{-50.5,-50}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(CO2_unbound.q_out, Hb0.CO2) annotation (Line(
            points={{68,40},{4,40},{4,82},{-4,82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb0.H, pH.q_out) annotation (Line(
            points={{-4,96.6},{10,96.6},{10,82},{24,82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb1.H, Hb0.H) annotation (Line(
            points={{-4,58.6},{10,58.6},{10,96.6},{-4,96.6}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb3.H, Hb0.H) annotation (Line(
            points={{-4,-25.4},{10,-25.4},{10,96.6},{-4,96.6}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb4.H, Hb0.H) annotation (Line(
            points={{-4,-69.4},{10,-69.4},{10,96.6},{-4,96.6}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb2.H, Hb0.H) annotation (Line(
            points={{-4,18.6},{10,18.6},{10,96.6},{-4,96.6}},
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
        connect(Hb0.Hbtn, K1.products[2]) annotation (Line(
            points={{-23.8,88},{-43.5,88},{-43.5,78}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb1.Hbtn, K1.substrates[1]) annotation (Line(
            points={{-23.8,50},{-44,50},{-44,58}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb1.Hbtn, K2.products[2]) annotation (Line(
            points={{-23.8,50},{-44,50},{-44,38},{-45.5,38}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb2.Hbtn, K2.substrates[1]) annotation (Line(
            points={{-23.8,10},{-46,10},{-46,18}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb2.Hbtn, K3.products[2]) annotation (Line(
            points={{-23.8,10},{-46,10},{-46,-8},{-47.5,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb3.Hbtn, K3.substrates[1]) annotation (Line(
            points={{-23.8,-34},{-48,-34},{-48,-28}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb3.Hbtn, K4.products[2]) annotation (Line(
            points={{-23.8,-34},{-48,-34},{-48,-50},{-49.5,-50}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb4.Hbtn, K4.substrates[1]) annotation (Line(
            points={{-23.8,-78},{-50,-78},{-50,-70}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb4.tHb_u, totalHemoglobin.fragment[1]) annotation (Line(
            points={{-24,-85.8},{-32,-85.8},{-32,-96},{22,-96},{22,10.4},{44,10.4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb3.tHb_u, totalHemoglobin.fragment[2]) annotation (Line(
            points={{-24,-41.8},{-32,-41.8},{-32,-48},{20,-48},{20,11.2},{44,11.2}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb2.tHb_u, totalHemoglobin.fragment[3]) annotation (Line(
            points={{-24,2.2},{-32,2.2},{-32,-4},{18,-4},{18,12},{44,12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb1.tHb_u, totalHemoglobin.fragment[4]) annotation (Line(
            points={{-24,42.2},{-28,42.2},{-28,34},{16,34},{16,12.8},{44,12.8}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb0.tHb_u, totalHemoglobin.fragment[5]) annotation (Line(
            points={{-24,80.2},{-28,80.2},{-28,64},{18,64},{18,13.6},{44,13.6}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(Hb1.tHb_u, sO2.u[4]) annotation (Line(
            points={{-24,42.2},{-28,42.2},{-28,34},{16,34},{16,-18.5},{60,-18.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb2.tHb_u, sO2.u[3]) annotation (Line(
            points={{-24,2.2},{-32,2.2},{-32,2},{-32,2},{-32,-4},{18,-4},{18,-19.5},{60,
                -19.5}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(Hb3.tHb_u, sO2.u[2]) annotation (Line(
            points={{-24,-41.8},{-32,-41.8},{-32,-48},{20,-48},{20,-20.5},{60,-20.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb4.tHb_u, sO2.u[1]) annotation (Line(
            points={{-24,-85.8},{-32,-85.8},{-32,-96},{22,-96},{22,-21.5},{60,-21.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(partialPressure1.q_out,oxygen_in_air. q_out)
                                                  annotation (Line(
            points={{-84,16},{-84,24}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(partialPressure1.q_in,oxygen_unbound. q_out) annotation (Line(
            points={{-84,-2},{-84,-18}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(Hb0.internalHeat, internalHeat.u[1]) annotation (Line(
            points={{-24,78.2},{-24,66},{34,66},{34,-53.6},{42,-53.6}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb1.internalHeat, internalHeat.u[2]) annotation (Line(
            points={{-24,40.2},{-26,40.2},{-26,36},{34,36},{34,-52.8},{42,-52.8}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb2.internalHeat, internalHeat.u[3]) annotation (Line(
            points={{-24,0.2},{-28,0.2},{-28,-2},{34,-2},{34,-52},{42,-52}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb3.internalHeat, internalHeat.u[4]) annotation (Line(
            points={{-24,-43.8},{-28,-43.8},{-28,-46},{34,-46},{34,-51.2},{42,
                -51.2}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(Hb4.internalHeat, internalHeat.u[5]) annotation (Line(
            points={{-24,-87.8},{-28,-87.8},{-28,-90},{34,-90},{34,-52},{42,-52},
                {42,-50.4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(gain.u, sO2.y) annotation (Line(
            points={{37.2,-88},{34,-88},{34,-96},{98,-96},{98,-20},{83,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(internalHeat.y, der1.u) annotation (Line(
            points={{65,-52},{68,-52},{68,-66},{48,-66},{48,-76},{51.2,-76}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(gain.y, der2.u) annotation (Line(
            points={{46.4,-88},{51.2,-88}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(clock.y, oxygen_in_air.partialPressure) annotation (Line(
            points={{-84,59},{-84,44}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(derInternalHeat_per_derO2.u2, der2.y) annotation (Line(
            points={{70,-88},{60.4,-88}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(derInternalHeat_per_derO2.u1, der1.y) annotation (Line(
            points={{70,-76},{60.4,-76}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(CO2_unbound.q_out, partialPressure2.q_in) annotation (Line(
            points={{68,40},{68,54}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2.q_out, partialPressure2.q_out) annotation (Line(
            points={{76,82},{68,82},{68,72}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        annotation (          experiment(
            StopTime=15000,
            Tolerance=1e-014,
            __Dymola_Algorithm="Euler"), Documentation(info="<html>
<p>Before silumation in &QUOT;Dymola 2014 FD01&QUOT; please set environment variable &QUOT;<code><b>Advanced.Define.NonLinearIterations&nbsp;=&nbsp;3&QUOT;</b></code> and chose &QUOT;Euler&QUOT; method!</p>

<p>[1] Mateják M, Kulhánek T, Matouaek S. Adair-Based Hemoglobin Equilibrium with Oxygen, Carbon Dioxide and Hydrogen Ion Activity. Scandinavian Journal of Clinical &AMP; Laboratory Investigation; 2015</p>

<p>[2] Bauer C, Schr&ouml;der E. Carbamino compounds of haemoglobin in human adult and foetal blood. The Journal of physiology 1972;227:457-71.</p>

<p>[3] Siggaard-Andersen O. Oxygen-Linked Hydrogen Ion Binding of Human Hemoglobin. Effects of Carbon Dioxide and 2, 3-Diphosphoglycerate I. Studies on Erythrolysate. Scandinavian Journal of Clinical &AMP; Laboratory Investigation 1971;27:351-60.</p>

<p>[4] Severinghaus JW. Simple, accurate equations for human blood O2 dissociation computations. Journal of Applied Physiology 1979;46:599-602.</p>
</html>", revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end Hemoglobin_MKM_Adair;

      package Develop
        extends Modelica.Icons.UnderConstruction;

        model QuaternaryForm
          "Model of hemoglobin space-structure form (can be parametrized as relaxed or tensed)"

          parameter Boolean isDependent = false;

          parameter Physiolibrary.Types.Concentration KA=10^(-6.89 + 3)
            "dissociation coefficient for acid chains of subunit";
          parameter Physiolibrary.Types.Concentration Kz=10^(-7.25 + 3)
            "valine 1 amino terminus dissociation coefficient of protonation to NH3+";
          parameter Physiolibrary.Types.Concentration Kc=10^(-8.35 + 3)
            "valine 1 amino terminus dissociation coefficient of protonation to NH3+";
          parameter Physiolibrary.Types.Concentration KO2=0.000671946
            "oxygen dissociation coefficient of hemoglobin subunit";

          Physiolibrary.Chemical.Components.Speciation Speciation(
              NumberOfSubunits=12)
            annotation (Placement(transformation(extent={{60,-20},{40,0}})));
          Physiolibrary.Chemical.Components.Substance OxyHm[4](
            each solute_start=0,
            each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            isDependent={isDependent,true,true,true})
            "Oxygenated subunit of hemoglobin tetramer" annotation (Placement(
                transformation(extent={{-90,-68},{-70,-48}})));
          Physiolibrary.Chemical.Components.ChemicalReaction oxygenation1[4](each
              nP=2, each K=KO2) annotation (Placement(transformation(extent={{-62,
                    -68},{-42,-48}})));
          Physiolibrary.Chemical.Components.Substance DeoxyHm[4](each
              Simulation=Physiolibrary.Types.SimulationType.SteadyState, each
              solute_start=1e-08) "Deoxygenated subunit of hemoglobin tetramer"
            annotation (Placement(transformation(extent={{-34,-68},{-14,-48}})));

          Modelica.Blocks.Math.Add add[4] annotation (Placement(transformation(
                extent={{-4,-4},{4,4}},
                rotation=270,
                origin={-58,-80})));
          Physiolibrary.Chemical.Interfaces.ChemicalPort_a O2 annotation (
              Placement(transformation(extent={{-26,-50},{-6,-30}}),
                iconTransformation(extent={{-26,-50},{-6,-30}})));
          Physiolibrary.Chemical.Interfaces.ChemicalPort_a sForm annotation (
              Placement(transformation(extent={{72,-54},{92,-34}}),
                iconTransformation(extent={{68,-50},{88,-30}})));
          Physiolibrary.Chemical.Interfaces.ChemicalPort_a H
            "hydrogen ion (proton)" annotation (Placement(transformation(extent=
                   {{-32,18},{-12,38}}), iconTransformation(extent={{-32,18},{-12,
                    38}})));
          Physiolibrary.Chemical.Components.Substance A[4](each Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, each
              solute_start=1e-08)
            "residual acid chains of hemoglobin subunits "
            annotation (Placement(transformation(extent={{-24,-14},{-4,6}})));
          Physiolibrary.Chemical.Components.Substance HA[4](
            each solute_start=0,
            each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            each isDependent=true)
            "residual acid chains of hemoglobin subunits "
            annotation (Placement(transformation(extent={{-90,-14},{-70,6}})));
          Physiolibrary.Chemical.Components.ChemicalReaction protonation1[4](each
              nP=2, each K=KA)
            annotation (Placement(transformation(extent={{-62,-14},{-42,6}})));
          Modelica.Blocks.Math.Add add1[
                                       4] annotation (Placement(transformation(
                extent={{-4,-4},{4,4}},
                rotation=270,
                origin={-52,-24})));
          Physiolibrary.Chemical.Components.Substance NH2[4](each Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, each
              solute_start=1e-08) "Val1 terminal of hemoglobin subunits "
            annotation (Placement(transformation(extent={{-10,52},{10,72}})));
          Physiolibrary.Chemical.Components.Substance NH3[4](
            each solute_start=0,
            each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            each isDependent=true) "Val1 terminal of hemoglobin subunits "
            annotation (Placement(transformation(extent={{-86,52},{-66,72}})));
          Physiolibrary.Chemical.Components.ChemicalReaction protonation2[4](each
              nP=2, each K=Kz)
            annotation (Placement(transformation(extent={{-58,52},{-38,72}})));
          Modelica.Blocks.Math.Add3 add2[
                                       4] annotation (Placement(transformation(
                extent={{-4,-4},{4,4}},
                rotation=270,
                origin={0,40})));
          Physiolibrary.Chemical.Interfaces.ChemicalPort_a CO2 annotation (
              Placement(transformation(extent={{-6,76},{14,96}}),
                iconTransformation(extent={{-6,76},{14,96}})));
          Physiolibrary.Chemical.Components.ChemicalReaction carboxylation[4](
            each nP=2,
            each nS=2,
            each K=Kc)
            "Carboxylation of Valin1 amino terminus of hemogloni subunit"
            annotation (Placement(transformation(extent={{36,52},{56,72}})));
          Physiolibrary.Chemical.Components.Substance NHCOO[4](each Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, each
              solute_start=1e-08) "Val1 terminal of hemoglobin subunits "
            annotation (Placement(transformation(extent={{66,52},{86,72}})));
          Physiolibrary.Types.RealIO.AmountOfSubstanceOutput tAmount(start=
                1e-08) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={50,-86}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={40,-90})));
          Physiolibrary.Types.RealIO.AmountOfSubstanceOutput protonation
            annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                  origin={100,14}), iconTransformation(extent={{-10,-10},{10,10}},
                  origin={90,20})));
          Modelica.Blocks.Math.Sum add3(k=cat(
                1,
                -ones(4),
                ones(8)), nin=12)
            annotation (Placement(transformation(extent={{78,10},{86,18}})));
          Modelica.Blocks.Math.Sum add4(nin=4)
            annotation (Placement(transformation(extent={{-4,-4},{4,4}},
                rotation=270,
                origin={-80,-86})));
          Physiolibrary.Types.RealIO.AmountOfSubstanceOutput oxygenation
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-80,-110}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-80,-90})));
        equation

          connect(OxyHm.solute, add.u2) annotation (Line(
              points={{-80,-68},{-80,-74},{-60.4,-74},{-60.4,-75.2}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(DeoxyHm.solute, add.u1) annotation (Line(
              points={{-24,-68},{-26,-68},{-26,-75.2},{-55.6,-75.2}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(OxyHm.q_out, oxygenation1.substrates[1]) annotation (Line(
              points={{-80,-58},{-62,-58}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(oxygenation1.products[1], DeoxyHm.q_out) annotation (Line(
              points={{-42,-58.5},{-34,-58.5},{-34,-58},{-24,-58}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(Speciation.specificForm, sForm) annotation (Line(
              points={{40,-18},{40,-44},{82,-44}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(oxygenation1[1].products[2], O2) annotation (Line(
              points={{-42,-57.5},{-36,-57.5},{-36,-40},{-16,-40}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(oxygenation1[2].products[2], O2) annotation (Line(
              points={{-42,-57.5},{-36,-57.5},{-36,-40},{-16,-40}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(oxygenation1[3].products[2], O2) annotation (Line(
              points={{-42,-57.5},{-36,-57.5},{-36,-40},{-16,-40}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(oxygenation1[4].products[2], O2) annotation (Line(
              points={{-42,-57.5},{-36,-57.5},{-36,-40},{-16,-40}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(HA.q_out, protonation1.substrates[1]) annotation (Line(
              points={{-80,-4},{-62,-4}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(protonation1.products[1], A.q_out) annotation (Line(
              points={{-42,-4.5},{-42,-4},{-14,-4}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(H, protonation1[1].products[2]) annotation (Line(
              points={{-22,28},{-32,28},{-32,-3.5},{-42,-3.5}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(protonation1[2].products[2], H) annotation (Line(
              points={{-42,-3.5},{-32,-3.5},{-32,28},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(protonation1[3].products[2], H) annotation (Line(
              points={{-42,-3.5},{-32,-3.5},{-32,28},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(protonation1[4].products[2], H) annotation (Line(
              points={{-42,-3.5},{-32,-3.5},{-32,28},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(HA.solute, add1.u2) annotation (Line(
              points={{-80,-14},{-80,-19.2},{-54.4,-19.2}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(A.solute, add1.u1) annotation (Line(
              points={{-14,-14},{-14,-16},{-49.6,-16},{-49.6,-19.2}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(DeoxyHm.q_out, Speciation.specificSubunitForm[1:4])
            annotation (Line(
              points={{-24,-58},{-24,-52},{28,-52},{28,0},{50,0},{50,-0.416667}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));

          connect(A.q_out, Speciation.specificSubunitForm[5:8]) annotation (
              Line(
              points={{-14,-4},{-14,0.25},{50,0.25}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(add.y, Speciation.amountOfSubunit[1:4]) annotation (Line(
              points={{-58,-84.4},{-58,-86},{14,-86},{14,-10.8333},{42,-10.8333}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(add1.y, Speciation.amountOfSubunit[5:8]) annotation (Line(
              points={{-52,-28.4},{-52,-28},{12,-28},{12,-9.5},{42,-9.5}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(NH3.q_out, protonation2.substrates[1]) annotation (Line(
              points={{-76,62},{-58,62}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(protonation2.products[1], NH2.q_out) annotation (Line(
              points={{-38,61.5},{0,61.5},{0,62}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(H, protonation2[1].products[2]) annotation (Line(
              points={{-22,28},{-22,62.5},{-38,62.5}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(protonation2[2].products[2], H) annotation (Line(
              points={{-38,62.5},{-22,62.5},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(protonation2[3].products[2], H) annotation (Line(
              points={{-38,62.5},{-22,62.5},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(protonation2[4].products[2], H) annotation (Line(
              points={{-38,62.5},{-22,62.5},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(add2.y, Speciation.amountOfSubunit[9:12]) annotation (Line(
              points={{0,35.6},{0,32},{12,32},{12,-8.16667},{42,-8.16667}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(NH2.q_out, Speciation.specificSubunitForm[9:12]) annotation (
              Line(
              points={{0,62},{0,58},{28,58},{28,0.916667},{50,0.916667}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(NH2.q_out, carboxylation.substrates[1]) annotation (Line(
              points={{0,62},{0,61.5},{36,61.5}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(carboxylation.products[1], NHCOO.q_out) annotation (Line(
              points={{56,61.5},{76,61.5},{76,62}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(carboxylation[1].products[2], H) annotation (Line(
              points={{56,62.5},{62,62.5},{62,28},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(carboxylation[1].substrates[2], CO2) annotation (Line(
              points={{36,62.5},{26,62.5},{26,86},{4,86}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(carboxylation[2].products[2], H) annotation (Line(
              points={{56,62.5},{62,62.5},{62,28},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(carboxylation[2].substrates[2], CO2) annotation (Line(
              points={{36,62.5},{26,62.5},{26,86},{4,86}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(carboxylation[3].products[2], H) annotation (Line(
              points={{56,62.5},{62,62.5},{62,28},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(carboxylation[3].substrates[2], CO2) annotation (Line(
              points={{36,62.5},{26,62.5},{26,86},{4,86}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(carboxylation[4].products[2], H) annotation (Line(
              points={{56,62.5},{62,62.5},{62,28},{-22,28}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(carboxylation[4].substrates[2], CO2) annotation (Line(
              points={{36,62.5},{26,62.5},{26,86},{4,86}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(NH3.solute, add2.u3) annotation (Line(
              points={{-76,52},{-76,44.8},{-3.2,44.8}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(NH2.solute, add2.u2) annotation (Line(
              points={{0,52},{0,44.8}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(NHCOO.solute, add2.u1) annotation (Line(
              points={{76,52},{76,44.8},{3.2,44.8}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(Speciation.amount, tAmount) annotation (Line(
              points={{50,-18},{50,-86}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(HA.solute, add3.u[9:12]) annotation (Line(
              points={{-80,-14},{-80,-20},{-98,-20},{-98,14.7333},{77.2,14.7333}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(NH3.solute, add3.u[5:8]) annotation (Line(
              points={{-76,52},{-76,14.2},{77.2,14.2}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(NHCOO.solute, add3.u[1:4]) annotation (Line(
              points={{76,52},{76,13.6667},{77.2,13.6667}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add3.y, protonation) annotation (Line(
              points={{86.4,14},{100,14}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(OxyHm.solute, add4.u) annotation (Line(
              points={{-80,-68},{-80,-81.2}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(add4.y, oxygenation) annotation (Line(
              points={{-80,-90.4},{-80,-110}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics={
                Ellipse(
                  extent={{-94,-44},{4,-80}},
                  lineColor={127,0,0},
                  fillColor={255,170,170},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-98,16},{8,-28}},
                  lineColor={0,0,255},
                  fillColor={179,254,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-94,80},{98,42}},
                  lineColor={127,127,0},
                  fillColor={213,255,170},
                  fillPattern=FillPattern.Solid)}),
            Documentation(revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>[1] Morrow J, Matthew J, Wittebort R, Gurd F. Carbon 13 resonances of 13CO2 carbamino adducts of alpha and beta chains in human adult hemoglobin. Journal of Biological Chemistry 1976;251:477-84.</p>
<p>[2] Bauer C, Schr&ouml;der E. Carbamino compounds of haemoglobin in human adult and foetal blood. The Journal of physiology 1972;227:457-71.</p>
<p>[3] Antonini E, Wyman J, Brunori M, Fronticelli C, Bucci E, Rossi-Fanelli A. Studies on the relations between molecular and functional properties of hemoglobin V. The influence of temperature on the Bohr effect in human and in horse hemoglobin. Journal of Biological Chemistry 1965;240:1096-103.</p>
</html>"),  Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}),
                graphics={
                Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={0,0,127},
                  fillColor={215,215,215},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-94,78},{98,40}},
                  lineColor={127,127,0},
                  fillColor={213,255,170},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-98,14},{8,-30}},
                  lineColor={0,0,255},
                  fillColor={179,254,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-94,-46},{4,-82}},
                  lineColor={127,0,0},
                  fillColor={255,170,170},
                  fillPattern=FillPattern.Solid),
                Text(
                  extent={{-240,-110},{40,-160}},
                  lineColor={0,0,0},
                  textString="%name")}));
        end QuaternaryForm;

        model Hemoglobin2 "Hemoglobin model"

         extends Physiolibrary.SteadyStates.Interfaces.SteadyStateSystem(
              Simulation=Physiolibrary.Types.SimulationType.SteadyState);

        //  parameter GasSolubility alpha =  0.0105 * 1e-3 "oxygen solubility in plasma";   // by Siggaard Andersen: 0.0105 (mmol/l)/kPa

          parameter Physiolibrary.Types.Fraction L=7.0529*10^6
            "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
          parameter Physiolibrary.Types.Fraction Ln=26884.8
            "quaternaly form ratio for specific stripped species of hemoglobin tetramer";
                                         //L*0.00381188                                                                     //"=L*(fnT/fnR)^4 for pH=7.2464 and CO2=0";
          parameter Physiolibrary.Types.Fraction c=0.00431555
            "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
          parameter Physiolibrary.Types.Concentration KR=0.000671946
            "oxygen dissociation on relaxed(R) hemoglobin subunit";
                                                                      //*7.875647668393782383419689119171e-5
                                                                    //10.500001495896 7.8756465463794e-05

          parameter Physiolibrary.Types.Concentration KT=KR/c
            "oxygen dissociation on tensed(T) hemoglobin subunit";

          parameter Physiolibrary.Types.AmountOfSubstance totalAmountOfHemoglobin=0.001;

          Physiolibrary.Chemical.Components.ChemicalReaction quaternaryForm(K=Ln)
            annotation (Placement(transformation(extent={{-16,26},{4,46}})));

          QuaternaryForm R(
            KO2=KR,
            KA=10^(-6.89 + 3),
            Kz=10^(-7.25 + 3),
            Kc=10^(-8.35 + 3),
            isDependent=true)
            annotation (Placement(transformation(extent={{-40,30},{-20,50}})));
          QuaternaryForm T(
            KO2=KT,
            KA=10^(-7.52 + 3),
            Kz=10^(-7.73 + 3),
            Kc=10^(-7.54 + 3))
            annotation (Placement(transformation(extent={{32,30},{12,50}})));
          Physiolibrary.SteadyStates.Components.MolarConservationLaw totalHb(
            n=2,
            Total=totalAmountOfHemoglobin,
            Simulation=Physiolibrary.Types.SimulationType.SteadyState)
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,2})));
          Physiolibrary.Chemical.Interfaces.ChemicalPort_a H "H+ (proton)"
            annotation (Placement(transformation(extent={{6,66},{26,86}}),
                iconTransformation(extent={{90,90},{110,110}})));
          Physiolibrary.Chemical.Interfaces.ChemicalPort_a CO2 "carbon dioxide"
            annotation (Placement(transformation(extent={{-22,54},{-2,74}}),
                iconTransformation(extent={{14,40},{34,60}})));
          Physiolibrary.Chemical.Interfaces.ChemicalPort_a O2 "oxygen"
            annotation (Placement(transformation(extent={{-54,78},{-34,98}}),
                iconTransformation(extent={{90,-10},{110,10}})));
          Physiolibrary.Types.RealIO.FractionOutput protonation annotation (
              Placement(transformation(extent={{-10,-10},{10,10}}, origin={100,
                    -40}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-80,-90})));
          Modelica.Blocks.Math.Add add(k1=1/4, k2=1/4)
            annotation (Placement(transformation(extent={{16,-40},{26,-30}})));
          Modelica.Blocks.Math.Division division
            annotation (Placement(transformation(extent={{42,-46},{54,-34}})));
          Modelica.Blocks.Math.Add add1(
                                       k1=1/4, k2=1/4)
            annotation (Placement(transformation(extent={{42,-62},{54,-50}})));
          Modelica.Blocks.Math.Division division1
            annotation (Placement(transformation(extent={{66,-86},{78,-74}})));
          Physiolibrary.Types.RealIO.FractionOutput sO2
            annotation (Placement(transformation(extent={{90,-90},{110,-70}})));
        equation

          connect(R.CO2, CO2) annotation (Line(
              points={{-29.6,48.6},{-29.6,64},{-12,64}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(CO2, T.CO2) annotation (Line(
              points={{-12,64},{22,64},{22,48},{21.6,48},{21.6,48.6}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));

          connect(R.O2, O2) annotation (Line(
              points={{-31.6,36},{-44,36},{-44,88}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(R.H, H) annotation (Line(
              points={{-32.2,42.8},{-32.2,42},{-34,42},{-34,76},{16,76}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(R.sForm, quaternaryForm.substrates[1]) annotation (Line(
              points={{-22.2,36},{-16,36}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(T.O2, O2) annotation (Line(
              points={{23.6,36},{36,36},{36,88},{-44,88}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(T.H, H) annotation (Line(
              points={{24.2,42.8},{26,42.8},{26,76},{16,76}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(T.sForm, quaternaryForm.products[1]) annotation (Line(
              points={{14.2,36},{4,36}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(R.tAmount, totalHb.fragment[1]) annotation (Line(
              points={{-26,31},{-26,18},{-4,18},{-4,12},{-5,12}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(T.tAmount, totalHb.fragment[2]) annotation (Line(
              points={{18,31},{18,18},{-3,18},{-3,12}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(T.protonation, add.u1) annotation (Line(
              points={{13,42},{13,-32},{15,-32}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(R.protonation, add.u2) annotation (Line(
              points={{-21,42},{-21,-38},{15,-38}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(division.u2, totalHb.totalAmountOfSubstance) annotation (Line(
              points={{40.8,-43.6},{-4,-43.6},{-4,-8}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(division.u1, add.y) annotation (Line(
              points={{40.8,-36.4},{38,-36.4},{38,-35},{26.5,-35}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(protonation, division.y) annotation (Line(
              points={{100,-40},{54.6,-40}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(division1.u1, add1.y) annotation (Line(
              points={{64.8,-76.4},{56,-76.4},{56,-56},{54.6,-56}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(R.oxygenation, add1.u2) annotation (Line(
              points={{-38,31},{-38,-59.6},{40.8,-59.6}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(T.oxygenation, add1.u1) annotation (Line(
              points={{30,31},{30,-52.4},{40.8,-52.4}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(totalHb.totalAmountOfSubstance, division1.u2) annotation (
              Line(
              points={{-4,-8},{-4,-83.6},{64.8,-83.6}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(division1.y, sO2) annotation (Line(
              points={{78.6,-80},{100,-80}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (             Documentation(revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>[1] Monod J, Wyman J, Changeux J-P. On the nature of allosteric transitions: a plausible model. Journal of Molecular Biology 1965;12:88-118.</p>
</html>"));
        end Hemoglobin2;

        model Hemoglobin_oxygenation "Hemoglobin oxygenation experiment"

          import Physiolibrary.Types.*;

         extends Modelica.Icons.Example;

          Physiolibrary.Chemical.Components.Substance oxygen_unbound(Simulation=
               SimulationType.SteadyState, solute_start=0.000001*
                7.875647668393782383419689119171e-5)
            annotation (Placement(transformation(extent={{-4,-2},{16,18}})));
          Modelica.Blocks.Sources.Clock clock(offset=1e-06)
            annotation (Placement(transformation(extent={{-40,74},{-20,94}})));
          Physiolibrary.Chemical.Sources.UnlimitedGasStorage oxygen_in_air(
            Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            usePartialPressureInput=true,
            T=310.15,
            isIsolatedInSteadyState=false) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={6,60})));
          Physiolibrary.Chemical.Components.GasSolubility partialPressure1(
            kH_T0(displayUnit="(mmol/l)/kPa at 25degC") = 0.026029047188736,
            T=310.15,
            C=1700) annotation (Placement(transformation(extent={{-10,-10},{10,
                    10}}, origin={6,32})));
          Physiolibrary.Chemical.Sources.UnlimitedSolutionStorage pH(
            q_out(conc(nominal=10^(-7.4 + 3))),
            Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            Conc=10^(-7.2464 + 3),
            isIsolatedInSteadyState=false) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=180,
                origin={62,-10})));
          Physiolibrary.Chemical.Sources.UnlimitedGasStorage CO2_gas(Simulation=
               Physiolibrary.Types.SimulationType.SteadyState, PartialPressure=
                5332.8954966) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-34,56})));
          Physiolibrary.Chemical.Components.GasSolubility gasSolubility(C=2400,
              kH_T0(displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)
            annotation (Placement(transformation(extent={{-44,20},{-24,40}})));
          Physiolibrary.Chemical.Components.Substance CO2_liquid(Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, isDependent=
                true)
            annotation (Placement(transformation(extent={{-44,-4},{-24,16}})));
          Hemoglobin2 hemoglobin
            annotation (Placement(transformation(extent={{-26,-74},{-6,-54}})));
        equation

          connect(partialPressure1.q_out, oxygen_in_air.q_out)
                                                    annotation (Line(
              points={{6,42},{6,50}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(partialPressure1.q_in, oxygen_unbound.q_out) annotation (Line(
              points={{6,24},{6,8}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(clock.y, oxygen_in_air.partialPressure) annotation (Line(
              points={{-19,84},{6,84},{6,70}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(gasSolubility.q_in,CO2_liquid. q_out) annotation (Line(
              points={{-34,22},{-34,6}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(CO2_gas.q_out,gasSolubility. q_out) annotation (Line(
              points={{-34,46},{-34,40}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(hemoglobin.CO2, CO2_liquid.q_out) annotation (Line(
              points={{-13.6,-59},{-13.6,-23.5},{-34,-23.5},{-34,6}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(hemoglobin.H, pH.q_out) annotation (Line(
              points={{-6,-54},{26,-54},{26,-10},{52,-10}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(hemoglobin.O2, oxygen_unbound.q_out) annotation (Line(
              points={{-6,-64},{-10,-64},{-10,8},{6,8}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          annotation (            experiment(
              StopTime=18000,
              Tolerance=1e-014,
              __Dymola_Algorithm="Euler"), Documentation(revisions=
                        "<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>[1] Severinghaus JW. Simple, accurate equations for human blood O2 dissociation computations. Journal of Applied Physiology 1979;46:599-602.</p>
</html>"));
        end Hemoglobin_oxygenation;

        model Hemoglobin_titration "Hemoglobin titration experiment"

          import Physiolibrary.Types.*;

         extends Modelica.Icons.Example;

         extends Physiolibrary.SteadyStates.Interfaces.SteadyStateSystem(
              Simulation=SimulationType.SteadyState);

        //  parameter GasSolubility alpha =  0.0105 * 1e-3 "oxygen solubility in plasma";   // by Siggaard Andersen: 0.0105 (mmol/l)/kPa

          parameter Fraction L = 7.0529*10^6
            "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
          parameter Fraction Ln = 26884.8
            "quaternaly form ratio for specific stripped species of hemoglobin tetramer";
                                         //L*0.00381188                                                                     //"=L*(fnT/fnR)^4 for pH=7.2464 and CO2=0";
          parameter Fraction c = 0.00431555
            "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
          parameter Concentration KR = 0.000671946
            "oxygen dissociation on relaxed(R) hemoglobin subunit";
                                                                      //*7.875647668393782383419689119171e-5
                                                                    //10.500001495896 7.8756465463794e-05

          parameter Concentration KT=KR/c
            "oxygen dissociation on tensed(T) hemoglobin subunit";

          parameter AmountOfSubstance totalAmountOfHemoglobin=0.001;

          Physiolibrary.Chemical.Components.Substance oxygen_unbound(Simulation=
               SimulationType.SteadyState, solute_start=0.000001*
                7.875647668393782383419689119171e-5)
            annotation (Placement(transformation(extent={{-4,-2},{16,18}})));
          Modelica.Blocks.Sources.Clock clock(offset=6.7)
            annotation (Placement(transformation(extent={{30,34},{50,54}})));
          Physiolibrary.Chemical.Sources.UnlimitedGasStorage oxygen_in_air(
            Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            usePartialPressureInput=false,
            PartialPressure=0,
            T=310.15,
            isIsolatedInSteadyState=false) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={6,60})));
          Physiolibrary.Chemical.Components.GasSolubility partialPressure1(
            kH_T0(displayUnit="(mmol/l)/kPa at 25degC") = 0.024913516594933,
            T=310.15,
            C=1700) annotation (Placement(transformation(extent={{-10,-10},{10,
                    10}}, origin={6,32})));
          Physiolibrary.Chemical.Sources.UnlimitedSolutionStorage pH(
            q_out(conc(nominal=10^(-7.4 + 3))),
            isIsolatedInSteadyState=false,
            Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            Conc=10^(-7.2464 + 3),
            useConcentrationInput=true) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=180,
                origin={62,-10})));
          Physiolibrary.Chemical.Sources.UnlimitedGasStorage CO2_gas(
            Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            usePartialPressureInput=false,
            PartialPressure=0) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-34,56})));
          Physiolibrary.Chemical.Components.GasSolubility gasSolubility(C=2400,
              kH_T0(displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)
            annotation (Placement(transformation(extent={{-44,20},{-24,40}})));
          Physiolibrary.Chemical.Components.Substance CO2_liquid(Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, isDependent=
                true)
            annotation (Placement(transformation(extent={{-44,-4},{-24,16}})));
          Hemoglobin2 deoxyhemoglobin
            annotation (Placement(transformation(extent={{-22,-68},{-2,-48}})));
          Physiolibrary.Types.RealIO.FractionOutput protonation
            "allosteric-dependent protonation"
            annotation (Placement(transformation(extent={{68,-76},{88,-56}})));
          Physiolibrary.Blocks.Math.Power pow annotation (Placement(
                transformation(
                extent={{-4,-4},{4,4}},
                rotation=270,
                origin={92,38})));
          Modelica.Blocks.Math.Gain gain(k=-1)
            annotation (Placement(transformation(extent={{62,34},{82,54}})));
          Modelica.Blocks.Math.Gain toMolPerM3(k=1000)
            "from mol/liter to mmol/liter (=mol/m3)" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={92,12})));
        equation

          connect(partialPressure1.q_out, oxygen_in_air.q_out)
                                                    annotation (Line(
              points={{6,42},{6,50}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(partialPressure1.q_in, oxygen_unbound.q_out) annotation (Line(
              points={{6,24},{6,8}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(gasSolubility.q_in,CO2_liquid. q_out) annotation (Line(
              points={{-34,22},{-34,6}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(CO2_gas.q_out,gasSolubility. q_out) annotation (Line(
              points={{-34,46},{-34,40}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(deoxyhemoglobin.CO2, CO2_liquid.q_out) annotation (Line(
              points={{-9.6,-53},{-9.6,-23.5},{-34,-23.5},{-34,6}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(clock.y, gain.u) annotation (Line(
              points={{51,44},{60,44}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(gain.y, pow.exponent) annotation (Line(
              points={{83,44},{90,44},{90,42},{89.6,42}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(pH.concentration, toMolPerM3.y) annotation (Line(
              points={{72,-10},{92,-10},{92,1}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(toMolPerM3.u, pow.y) annotation (Line(
              points={{92,24},{92,33.6}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(deoxyhemoglobin.H, pH.q_out) annotation (Line(
              points={{-2,-48},{26,-48},{26,-10},{52,-10}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(deoxyhemoglobin.O2, oxygen_unbound.q_out) annotation (Line(
              points={{-2,-58},{-6,-58},{-6,8},{6,8}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(deoxyhemoglobin.protonation, protonation) annotation (Line(
              points={{-20,-67},{-20,-66},{78,-66}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (            experiment(StopTime=1.3), Documentation(revisions=
                        "<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>[1] Siggaard-Andersen O, Garby L. The Bohr effect and the Haldane effect. Scandinavian Journal of Clinical &AMP; Laboratory Investigation 1973;31:1-8.</p>
</html>"));
        end Hemoglobin_titration;

        model Hemoglobin_titration_shifts
          "Hemoglobin titration shift caused by full deoxygenation (Bohr protons binding)"
          extends Modelica.Icons.Example;

          Hemoglobin_titration hemoglobin_titration
            annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
          Hemoglobin_titration hemoglobin_titration1(CO2_gas(PartialPressure(
                  displayUnit="kPa") = 1470))
            annotation (Placement(transformation(extent={{-28,60},{-8,80}})));
          Hemoglobin_titration hemoglobin_titration2(CO2_gas(PartialPressure(
                  displayUnit="kPa") = 4530))
            annotation (Placement(transformation(extent={{0,60},{20,80}})));
          Hemoglobin_titration hemoglobin_titration3(CO2_gas(PartialPressure(
                  displayUnit="kPa") = 10670))
            annotation (Placement(transformation(extent={{30,60},{50,80}})));
          Hemoglobin_titration hemoglobin_titration4(CO2_gas(PartialPressure(
                  displayUnit="kPa") = 26660))
            annotation (Placement(transformation(extent={{60,60},{80,80}})));
          Hemoglobin_titration hemoglobin_titration5(oxygen_in_air(
                PartialPressure=19998.35811225))
            annotation (Placement(transformation(extent={{-60,-26},{-40,-6}})));
          Hemoglobin_titration hemoglobin_titration6(oxygen_in_air(
                PartialPressure=19998.35811225), CO2_gas(PartialPressure(
                  displayUnit="kPa") = 1470))
            annotation (Placement(transformation(extent={{-28,-26},{-8,-6}})));
          Hemoglobin_titration hemoglobin_titration7(oxygen_in_air(
                PartialPressure=19998.35811225), CO2_gas(PartialPressure(
                  displayUnit="kPa") = 4530))
            annotation (Placement(transformation(extent={{0,-26},{20,-6}})));
          Hemoglobin_titration hemoglobin_titration8(oxygen_in_air(
                PartialPressure=19998.35811225), CO2_gas(PartialPressure(
                  displayUnit="kPa") = 10670))
            annotation (Placement(transformation(extent={{30,-26},{50,-6}})));
          Hemoglobin_titration hemoglobin_titration9(oxygen_in_air(
                PartialPressure=19998.35811225), CO2_gas(PartialPressure(
                  displayUnit="kPa") = 26660))
            annotation (Placement(transformation(extent={{60,-26},{80,-6}})));
          Modelica.Blocks.Math.Feedback dH
            annotation (Placement(transformation(extent={{-54,22},{-34,42}})));
          Modelica.Blocks.Math.Feedback dH1
            annotation (Placement(transformation(extent={{-26,22},{-6,42}})));
          Modelica.Blocks.Math.Feedback dH2
            annotation (Placement(transformation(extent={{10,22},{30,42}})));
          Modelica.Blocks.Math.Feedback dH3
            annotation (Placement(transformation(extent={{36,22},{56,42}})));
          Modelica.Blocks.Math.Feedback dH4
            annotation (Placement(transformation(extent={{70,20},{90,40}})));
        equation
          connect(hemoglobin_titration.protonation, dH.u1) annotation (Line(
              points={{-42.2,63.4},{-42.2,47.7},{-52,47.7},{-52,32}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(dH.u2, hemoglobin_titration5.protonation) annotation (Line(
              points={{-44,24},{-42,24},{-42,-22.6},{-42.2,-22.6}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(hemoglobin_titration1.protonation, dH1.u1) annotation (Line(
              points={{-10.2,63.4},{-10.2,47.7},{-24,47.7},{-24,32}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(dH1.u2, hemoglobin_titration6.protonation) annotation (Line(
              points={{-16,24},{-14,24},{-14,-22.6},{-10.2,-22.6}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(dH2.u2, hemoglobin_titration7.protonation) annotation (Line(
              points={{20,24},{20,-22.6},{17.8,-22.6}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(hemoglobin_titration2.protonation, dH2.u1) annotation (Line(
              points={{17.8,63.4},{17.8,47.7},{12,47.7},{12,32}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(hemoglobin_titration3.protonation, dH3.u1) annotation (Line(
              points={{47.8,63.4},{47.8,46.7},{38,46.7},{38,32}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(dH3.u2, hemoglobin_titration8.protonation) annotation (Line(
              points={{46,24},{47.8,24},{47.8,-22.6}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(hemoglobin_titration4.protonation, dH4.u1) annotation (Line(
              points={{77.8,63.4},{77.8,46.7},{72,46.7},{72,30}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(dH4.u2, hemoglobin_titration9.protonation) annotation (Line(
              points={{80,22},{80,-22.6},{77.8,-22.6}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (
            experiment(
              StopTime=1.1,
              Tolerance=1e-014,
              __Dymola_Algorithm="Euler"), Documentation(info=
                   "<html>
<p>[1] Siggaard-Andersen O. Oxygen-Linked Hydrogen Ion Binding of Human Hemoglobin. Effects of Carbon Dioxide and 2, 3-Diphosphoglycerate I. Studies on Erythrolysate. Scandinavian Journal of Clinical &AMP; Laboratory Investigation 1971;27:351-60.</p>
</html>"));
        end Hemoglobin_titration_shifts;
      end Develop;

      model Allosteric_Hemoglobin_MWC2 "Monod,Wyman,Changeux (1965)"
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
        constant Real KRx=KR*AmountOfSolutionIn1L "Mole fraction based KR";

      //Relative Gibbs formation energies of the substances in the system:
        constant Real RT = Modelica.Constants.R*T;
        constant Modelica.SIunits.MolarEnergy
          GO2aq=-RT*log(0.0013*0.018),
          GR0=0,                            GT0=GR0 -RT*log(L),
          GR1=GR0+GO2aq +RT*log(KRx/4),     GT1=GR1 -RT*log(c*L),
          GR2=GR1+GO2aq +RT*log(KRx/(3/2)), GT2=GR2 -RT*log(c^2*L),
          GR3=GR2+GO2aq +RT*log(KRx/(2/3)), GT3=GR3 -RT*log(c^3*L),
          GR4=GR3+GO2aq +RT*log(KRx*4),     GT4=GR4 -RT*log(c^4*L);

        Sources.AmbientMolality
                             oxygen_unbound(substance( DfG_25degC=GO2aq),
          Molality=0.0001,
          AmountOfSolutionPer1kgOfSolvent=AmountOfSolutionIn1L)
          annotation (Placement(transformation(extent={{-56,-44},{-36,-24}})));

        Components.Substance T0(substance( DfG_25degC=GT0), amountOfSubstance_start=
              THb)
          annotation (Placement(transformation(extent={{34,78},{54,98}})));

        Components.Substance T1(substance( DfG_25degC=GT1),
            amountOfSubstance_start=1e-10)
          annotation (Placement(transformation(extent={{34,36},{54,56}})));

        Components.Substance T2(substance( DfG_25degC=GT2),
            amountOfSubstance_start=1e-18)
          annotation (Placement(transformation(extent={{34,-10},{54,10}})));

        Components.Substance R1(substance( DfG_25degC=GR1),
            amountOfSubstance_start=1e-15)
          annotation (Placement(transformation(extent={{-20,36},{0,56}})));

        Components.Substance R2(substance( DfG_25degC=GR2),
            amountOfSubstance_start=1e-20)
          annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

        Components.Substance T3(substance( DfG_25degC=GT3),
            amountOfSubstance_start=1e-25)
          annotation (Placement(transformation(extent={{34,-54},{54,-34}})));

        Components.Substance R3(substance( DfG_25degC=GR3),
            amountOfSubstance_start=1e-25)
          annotation (Placement(transformation(extent={{-20,-54},{0,-34}})));

        Components.Substance T4(substance( DfG_25degC=GT4),
            amountOfSubstance_start=1e-33)
          annotation (Placement(transformation(extent={{34,-92},{54,-72}})));

        Components.Substance R4(substance( DfG_25degC=GR4),
            amountOfSubstance_start=1e-31)
          annotation (Placement(transformation(extent={{-20,-92},{0,-72}})));

        Components.Substance R0(substance( DfG_25degC=GR0),
            amountOfSubstance_start=1e-10)
          annotation (Placement(transformation(extent={{-20,78},{0,98}})));

        Components.Reaction quaternaryForm annotation (Placement(transformation(extent={{4,78},{24,98}})));
        Components.Reaction oxyR1(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,64})));
        Components.Reaction oxyT1(nP=2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,64})));
        Components.Reaction oxyR2(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,22})));
        Components.Reaction oxyR3(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-24})));
        Components.Reaction oxyR4(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-66})));
        Components.Reaction oxyT2(nP=2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,22})));
        Components.Reaction oxyT3(nP=2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-24})));
        Components.Reaction oxyT4(nP=2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-66})));

        Modelica.Blocks.Math.Sum oxygen_bound(k={1,1,2,2,3,3,4,4}, nin=8)
          annotation (Placement(transformation(extent={{72,-42},{82,-32}})));
        Modelica.Blocks.Math.Division sO2_ "hemoglobin oxygen saturation"
          annotation (Placement(transformation(extent={{86,-60},{96,-50}})));
        Modelica.Blocks.Math.Sum tHb(nin=10, k=4*ones(10))
          annotation (Placement(transformation(extent={{70,-80},{80,-70}})));

        Components.Solution solution(amountOfSolution_start=
              AmountOfSolutionIn1L, Isothermal=false)
          annotation (Placement(transformation(extent={{-56,-102},{100,104}})));

      equation
       //  sO2 = (R1.amountOfSubstance + 2*R2.amountOfSubstance + 3*R3.amountOfSubstance + 4*R4.amountOfSubstance + T1.amountOfSubstance + 2*T2.amountOfSubstance + 3*T3.amountOfSubstance + 4*T4.amountOfSubstance)/(4*totalAmountOfHemoglobin);
      //   totalAmountOfRforms = R0.amountOfSubstance + R1.amountOfSubstance + R2.amountOfSubstance + R3.amountOfSubstance + R4.amountOfSubstance;
      //   totalAmountOfTforms = T0.amountOfSubstance + T1.amountOfSubstance + T2.amountOfSubstance + T3.amountOfSubstance + T4.amountOfSubstance;

      //   totalAmountOfHemoglobin*normalizedState[1] = totalAmountOfRforms + totalAmountOfTforms;

        connect(quaternaryForm.products[1],T0. port_a) annotation (Line(
            points={{24,88},{44,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR1.substrates[1],R1. port_a) annotation (Line(
            points={{-10,54},{-10,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R1.port_a,oxyR2. products[1]) annotation (Line(
            points={{-10,46},{-10,32},{-10.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR2.substrates[1],R2. port_a) annotation (Line(
            points={{-10,12},{-10,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.substrates[1],R3. port_a) annotation (Line(
            points={{-10,-34},{-10,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.products[1],R2. port_a) annotation (Line(
            points={{-10.5,-14},{-10.5,-7},{-10,-7},{-10,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R3.port_a,oxyR4. products[1]) annotation (Line(
            points={{-10,-44},{-10,-56},{-10.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR4.substrates[1],R4. port_a) annotation (Line(
            points={{-10,-76},{-10,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT1.products[1],T0. port_a) annotation (Line(
            points={{44.5,74},{44.5,88},{44,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT1.substrates[1],T1. port_a) annotation (Line(
            points={{44,54},{44,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(T1.port_a,oxyT2. products[1]) annotation (Line(
            points={{44,46},{44,32},{44.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT3.substrates[1],T3. port_a) annotation (Line(
            points={{44,-34},{44,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(T3.port_a,oxyT4. products[1]) annotation (Line(
            points={{44,-44},{44,-56},{44.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT4.substrates[1],T4. port_a) annotation (Line(
            points={{44,-76},{44,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,quaternaryForm. substrates[1]) annotation (Line(
            points={{-10,88},{4,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,oxyR1. products[1]) annotation (Line(
            points={{-10,88},{-10,74},{-10.5,74}},
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
        connect(oxyT2.substrates[1], T2.port_a) annotation (Line(
            points={{44,12},{44,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(T2.port_a, oxyT3.products[1]) annotation (Line(
            points={{44,0},{44,-14},{44.5,-14}},
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
        annotation (          experiment(
            StopTime=1e-005,
            Tolerance=1e-005,
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
<p><i>2013</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics),
          __Dymola_experimentSetupOutput);
      end Allosteric_Hemoglobin_MWC2;

      model Allosteric_Hemoglobin_MWC3 "Monod,Wyman,Changeux (1965)"
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
        constant Real KRx=KR*AmountOfSolutionIn1L "Mole fraction based KR";

      //Relative Gibbs formation energies of the substances in the system:
        constant Real RT = Modelica.Constants.R*T;
        constant Modelica.SIunits.MolarEnergy
          GO2aq=-RT*log(0.0013*0.018),
          GR0=0,                            GT0=GR0 -RT*log(L),
          GR1=GR0+GO2aq +RT*log(KRx/4),     GT1=GR1 -RT*log(c*L),
          GR2=GR1+GO2aq +RT*log(KRx/(3/2)), GT2=GR2 -RT*log(c^2*L),
          GR3=GR2+GO2aq +RT*log(KRx/(2/3)), GT3=GR3 -RT*log(c^3*L),
          GR4=GR3+GO2aq +RT*log(KRx*4),     GT4=GR4 -RT*log(c^4*L);

        Sources.AmbientMolality
                             oxygen_unbound(substance( DfG_25degC=GO2aq),
          AmountOfSolutionPer1kgOfSolvent=AmountOfSolutionIn1L,
          Molality=1e-6)
          annotation (Placement(transformation(extent={{-56,-44},{-36,-24}})));

        Components.Substance T0(substance( DfG_25degC=GT0), amountOfSubstance_start=
              THb)
          annotation (Placement(transformation(extent={{34,78},{54,98}})));

        Components.Substance T1(substance( DfG_25degC=GT1),
            amountOfSubstance_start=1e-10)
          annotation (Placement(transformation(extent={{34,36},{54,56}})));

        Components.Substance T2(substance( DfG_25degC=GT2),
            amountOfSubstance_start=1e-18)
          annotation (Placement(transformation(extent={{34,-10},{54,10}})));

        Components.Substance R1(substance( DfG_25degC=GR1),
            amountOfSubstance_start=1e-15)
          annotation (Placement(transformation(extent={{-20,36},{0,56}})));

        Components.Substance R2(substance( DfG_25degC=GR2),
            amountOfSubstance_start=1e-20)
          annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

        Components.Substance T3(substance( DfG_25degC=GT3),
            amountOfSubstance_start=1e-25)
          annotation (Placement(transformation(extent={{34,-54},{54,-34}})));

        Components.Substance R3(substance( DfG_25degC=GR3),
            amountOfSubstance_start=1e-25)
          annotation (Placement(transformation(extent={{-20,-54},{0,-34}})));

        Components.Substance T4(substance( DfG_25degC=GT4),
            amountOfSubstance_start=1e-33)
          annotation (Placement(transformation(extent={{34,-92},{54,-72}})));

        Components.Substance R4(substance( DfG_25degC=GR4),
            amountOfSubstance_start=1e-31)
          annotation (Placement(transformation(extent={{-20,-92},{0,-72}})));

        Components.Substance R0(substance( DfG_25degC=GR0),
            amountOfSubstance_start=1e-10)
          annotation (Placement(transformation(extent={{-20,78},{0,98}})));

        Components.Reaction quaternaryForm annotation (Placement(transformation(extent={{4,78},{24,98}})));
        Components.Reaction oxyR1(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,64})));
        Components.Reaction oxyR2(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,22})));
        Components.Reaction oxyR3(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-24})));
        Components.Reaction oxyR4(nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-66})));

        Modelica.Blocks.Math.Sum oxygen_bound(k={1,1,2,2,3,3,4,4}, nin=8)
          annotation (Placement(transformation(extent={{72,-42},{82,-32}})));
        Modelica.Blocks.Math.Division sO2_ "hemoglobin oxygen saturation"
          annotation (Placement(transformation(extent={{86,-60},{96,-50}})));
        Modelica.Blocks.Math.Sum tHb(nin=10, k=4*ones(10))
          annotation (Placement(transformation(extent={{70,-80},{80,-70}})));

        Components.Solution solution(amountOfSolution_start=
              AmountOfSolutionIn1L, Isothermal=false)
          annotation (Placement(transformation(extent={{-56,-102},{100,104}})));

        Components.Reaction quaternaryForm1 annotation (Placement(transformation(extent={{8,36},{
                  28,56}})));
        Components.Reaction quaternaryForm2 annotation (Placement(transformation(extent={{8,-10},
                  {28,10}})));
        Components.Reaction quaternaryForm3 annotation (Placement(transformation(extent={{8,-54},
                  {28,-34}})));
        Components.Reaction quaternaryForm4 annotation (Placement(transformation(extent={{10,-92},
                  {30,-72}})));
      equation
       //  sO2 = (R1.amountOfSubstance + 2*R2.amountOfSubstance + 3*R3.amountOfSubstance + 4*R4.amountOfSubstance + T1.amountOfSubstance + 2*T2.amountOfSubstance + 3*T3.amountOfSubstance + 4*T4.amountOfSubstance)/(4*totalAmountOfHemoglobin);
      //   totalAmountOfRforms = R0.amountOfSubstance + R1.amountOfSubstance + R2.amountOfSubstance + R3.amountOfSubstance + R4.amountOfSubstance;
      //   totalAmountOfTforms = T0.amountOfSubstance + T1.amountOfSubstance + T2.amountOfSubstance + T3.amountOfSubstance + T4.amountOfSubstance;

      //   totalAmountOfHemoglobin*normalizedState[1] = totalAmountOfRforms + totalAmountOfTforms;

        connect(quaternaryForm.products[1],T0. port_a) annotation (Line(
            points={{24,88},{44,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR1.substrates[1],R1. port_a) annotation (Line(
            points={{-10,54},{-10,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R1.port_a,oxyR2. products[1]) annotation (Line(
            points={{-10,46},{-10,32},{-10.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR2.substrates[1],R2. port_a) annotation (Line(
            points={{-10,12},{-10,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.substrates[1],R3. port_a) annotation (Line(
            points={{-10,-34},{-10,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.products[1],R2. port_a) annotation (Line(
            points={{-10.5,-14},{-10.5,-7},{-10,-7},{-10,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R3.port_a,oxyR4. products[1]) annotation (Line(
            points={{-10,-44},{-10,-56},{-10.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR4.substrates[1],R4. port_a) annotation (Line(
            points={{-10,-76},{-10,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,quaternaryForm. substrates[1]) annotation (Line(
            points={{-10,88},{4,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,oxyR1. products[1]) annotation (Line(
            points={{-10,88},{-10,74},{-10.5,74}},
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
        connect(R1.port_a,quaternaryForm1. substrates[1]) annotation (Line(
            points={{-10,46},{8,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm1.products[1],T1. port_a) annotation (Line(
            points={{28,46},{44,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R2.port_a,quaternaryForm2. substrates[1]) annotation (Line(
            points={{-10,0},{8,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm2.products[1], T2.port_a) annotation (Line(
            points={{28,0},{44,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(R3.port_a,quaternaryForm3. substrates[1]) annotation (Line(
            points={{-10,-44},{8,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm3.products[1],T3. port_a) annotation (Line(
            points={{28,-44},{44,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R4.port_a,quaternaryForm4. substrates[1]) annotation (Line(
            points={{-10,-82},{10,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm4.products[1],T4. port_a) annotation (Line(
            points={{30,-82},{44,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        annotation (          experiment(
            StopTime=1e-005,
            Tolerance=1e-005,
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
<p><i>2013</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics),
          __Dymola_experimentSetupOutput);
      end Allosteric_Hemoglobin_MWC3;
    end Hemoglobin;

    package AcidBase
      model WaterSelfIonization
        "H2O  <->  OH-   +   H+ (It is better to solve this model using Euler solver, because there is only time dependence/no integration needed/)"
        import Chemical;
          extends Modelica.Icons.Example;
        Components.Substance H3O(
          substance=Chemical.Examples.Substances.Hydronium_aqueous,
          AmountOfSolution(displayUnit="mol") = 1,
          amountOfSubstance_start=1e-7,
          useSolutionPort=false) annotation (Placement(transformation(extent={{
                  -10,-10},{10,10}}, origin={36,72})));
        Components.Substance OH(
          substance=Chemical.Examples.Substances.Hydroxide_aqueous,
          AmountOfSolution(displayUnit="mol") = 1,
          amountOfSubstance_start=1e-7,
          useSolutionPort=false) annotation (Placement(transformation(extent={{
                  -10,-10},{10,10}}, origin={36,26})));
        Components.Substance H2O(
          substance=Chemical.Examples.Substances.Water_liquid,
          AmountOfSolution(displayUnit="mol") = 1,
          amountOfSubstance_start=1,
          useSolutionPort=false) annotation (Placement(transformation(extent={{
                  -10,-10},{10,10}}, origin={-38,46})));
        Chemical.Components.Reaction waterDissociation(nP=2, s={2})
          annotation (Placement(transformation(extent={{-12,36},{8,56}})));
              Real pH, pH_;
        Components.Substance H_(
          AmountOfSolution(displayUnit="mol") = 1,
          amountOfSubstance_start=1e-7,
          substance=Chemical.Examples.Substances.Proton_aqueous,
          useSolutionPort=false) annotation (Placement(transformation(extent={{
                  -10,-10},{10,10}}, origin={34,-30})));
        Components.Substance OH_(
          substance=Chemical.Examples.Substances.Hydroxide_aqueous,
          AmountOfSolution(displayUnit="mol") = 1,
          amountOfSubstance_start=1e-7,
          useSolutionPort=false) annotation (Placement(transformation(extent={{
                  -10,-10},{10,10}}, origin={34,-76})));
        Components.Substance H2O_(
          substance=Chemical.Examples.Substances.Water_liquid,
          AmountOfSolution(displayUnit="mol") = 1,
          amountOfSubstance_start=1,
          useSolutionPort=false) annotation (Placement(transformation(extent={{
                  -10,-10},{10,10}}, origin={-40,-56})));
        Chemical.Components.Reaction waterDissociation_(nP=2)
          annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));
      equation
        pH = -log10(H3O.port_a.activityCoefficient*H3O.amountOfSubstance);
        pH_ = -log10(H_.port_a.activityCoefficient*H_.amountOfSubstance);

        connect(OH.port_a, waterDissociation.products[1]) annotation (Line(
            points={{36,26},{22,26},{22,45.5},{8,45.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(waterDissociation.products[2], H3O.port_a) annotation (Line(
            points={{8,46.5},{22,46.5},{22,72},{36,72}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.port_a, waterDissociation.substrates[1]) annotation (Line(
            points={{-38,46},{-12,46}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(OH_.port_a,waterDissociation_. products[1]) annotation (Line(
            points={{34,-76},{20,-76},{20,-56.5},{6,-56.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(waterDissociation_.products[2], H_.port_a) annotation (Line(
            points={{6,-55.5},{20,-55.5},{20,-30},{34,-30}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O_.port_a,waterDissociation_. substrates[1]) annotation (Line(
            points={{-40,-56},{-14,-56}},
            color={158,66,200},
            thickness=1,
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

      model WaterSelfIonization1
        "H2O  <->  OH-   +   H+ (It is better to solve this model using Euler solver, because there is only time dependence/no integration needed/)"
        import Chemical;
          extends Modelica.Icons.Example;
        Components.Substance H3O(
          substance=Chemical.Examples.Substances.Hydronium_aqueous,
          amountOfSubstance_start=1e-7,
          useSolutionPort=true) annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={36,72})));
        Components.Substance OH(
          substance=Chemical.Examples.Substances.Hydroxide_aqueous,
          amountOfSubstance_start=1e-7,
          useSolutionPort=true) annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={36,26})));
        Components.Substance H2O(
          substance=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=1,
          useSolutionPort=true) annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={-38,46})));
        Chemical.Components.Reaction waterDissociation(nP=2, s={2})
          annotation (Placement(transformation(extent={{-12,36},{8,56}})));
              Real pH, pH_;
        Components.Substance H_(
          amountOfSubstance_start=1e-7,
          substance=Chemical.Examples.Substances.Proton_aqueous,
          useSolutionPort=true) annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={34,-30})));
        Components.Substance OH_(
          substance=Chemical.Examples.Substances.Hydroxide_aqueous,
          amountOfSubstance_start=1e-7,
          useSolutionPort=true) annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={34,-76})));
        Components.Substance H2O_(
          substance=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=1,
          useSolutionPort=true) annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={-40,-56})));
        Chemical.Components.Reaction waterDissociation_(nP=2)
          annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));
        Chemical.Components.Solution solution(amountOfSolution_start=1)
          annotation (Placement(transformation(extent={{-74,4},{74,98}})));
        Chemical.Components.Solution solution1(amountOfSolution_start=1)
          annotation (Placement(transformation(extent={{-76,-94},{72,0}})));
      equation
        pH = -log10(H3O.port_a.activityCoefficient*H3O.amountOfSubstance);
        pH_ = -log10(H_.port_a.activityCoefficient*H_.amountOfSubstance);

        connect(OH.port_a, waterDissociation.products[1]) annotation (Line(
            points={{36,26},{22,26},{22,45.5},{8,45.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(waterDissociation.products[2], H3O.port_a) annotation (Line(
            points={{8,46.5},{22,46.5},{22,72},{36,72}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.port_a, waterDissociation.substrates[1]) annotation (Line(
            points={{-38,46},{-12,46}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(OH_.port_a,waterDissociation_. products[1]) annotation (Line(
            points={{34,-76},{20,-76},{20,-56.5},{6,-56.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(waterDissociation_.products[2], H_.port_a) annotation (Line(
            points={{6,-55.5},{20,-55.5},{20,-30},{34,-30}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O_.port_a,waterDissociation_. substrates[1]) annotation (Line(
            points={{-40,-56},{-14,-56}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.solutionPort, solution.port_a) annotation (Line(
            points={{-44,36},{-59.2,36},{-59.2,4}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(OH.solutionPort, solution.port_a) annotation (Line(
            points={{30,16},{30,4},{-59.2,4}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H3O.solutionPort, solution.port_a) annotation (Line(
            points={{30,62},{52,62},{52,4},{-59.2,4}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H2O_.solutionPort, solution1.port_a) annotation (Line(
            points={{-46,-66},{-61.2,-66},{-61.2,-94}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(OH_.solutionPort, solution1.port_a) annotation (Line(
            points={{28,-86},{28,-94},{-61.2,-94}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H_.solutionPort, solution1.port_a) annotation (Line(
            points={{28,-40},{52,-40},{52,-94},{-61.2,-94}},
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
      end WaterSelfIonization1;

      model CarbonDioxideInWater "CO2 as alone acid-base buffer"
        import Chemical;
          extends Modelica.Icons.Example;
        Components.Substance HCO3(substance=Chemical.Examples.Substances.Bicarbonate_aqueous)
          annotation (Placement(transformation(extent={{-18,46},{2,66}})));
        Chemical.Components.Reaction HendersonHasselbalch(nP=2, nS=2,
          Tau=1200) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-58,22},{-38,42}})));
        Sources.AirSubstance CO2_gas(
          substance=Chemical.Examples.Substances.CarbonDioxide_gas,
          PartialPressure=5332.8954966,
          TotalPressure=101325.0144354) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-80,84})));
        Components.Substance H(substance=Chemical.Examples.Substances.Proton_aqueous,
            amountOfSubstance_start=3e-8) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}}, origin={-8,12})));
        Components.GasSolubility gasSolubility(Tau=1200)
          annotation (Placement(transformation(extent={{-90,46},{-70,66}})));
                                              /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
        Components.Substance CO2_liquid(substance=Chemical.Examples.Substances.CarbonDioxide_aqueous)
          annotation (Placement(transformation(extent={{-90,22},{-70,42}})));
        Components.Substance CO3(substance=Chemical.Examples.Substances.Carbonate_aqueous)
          annotation (Placement(transformation(extent={{64,54},{84,74}})));
        Chemical.Components.Reaction c2(nP=2, nS=1)
          "K=10^(-10.33 + 3), dH=14.9kJ/mol"
          annotation (Placement(transformation(extent={{16,46},{36,66}})));
        Chemical.Components.Substance H2O(substance=Chemical.Examples.Substances.Water_liquid,
            amountOfSubstance_start=55.507)
          annotation (Placement(transformation(extent={{-90,-22},{-70,-2}})));
        Real pH;
      equation
        pH = -log10(H.port_a.activityCoefficient*H.amountOfSubstance);

        connect(CO2_gas.port_a, gasSolubility.gas_port) annotation (Line(
            points={{-80,74},{-80,66}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(gasSolubility.liquid_port, CO2_liquid.port_a) annotation (Line(
            points={{-80,46},{-80,32}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch.products[1], H.port_a) annotation (Line(
            points={{-38,31.5},{-24,31.5},{-24,12},{-8,12}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch.products[2], HCO3.port_a) annotation (Line(
            points={{-38,32.5},{-24,32.5},{-24,56},{-8,56}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HCO3.port_a, c2.substrates[1]) annotation (Line(
            points={{-8,56},{16,56}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(c2.products[1], H.port_a) annotation (Line(
            points={{36,55.5},{48,55.5},{48,12},{-8,12}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(c2.products[2], CO3.port_a) annotation (Line(
            points={{36,56.5},{52,56.5},{52,64},{74,64}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.port_a, HendersonHasselbalch.substrates[2]) annotation (
            Line(
            points={{-80,32},{-70,32},{-70,32.5},{-58,32.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.port_a, HendersonHasselbalch.substrates[1]) annotation (Line(
            points={{-80,-12},{-64,-12},{-64,31.5},{-58,31.5}},
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
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isDependent=true,
          solute_start=55.6*10^(-7.4),
          q_out(conc(nominal=10^(-7.4))),
          substance=Chemical.Examples.Substances.Proton)
          "hydrogen ions activity" annotation (Placement(transformation(extent=
                  {{-10,-10},{10,10}}, origin={36,-12})));

        Components.Substance H3PO4(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isDependent=true,
          solute_start=1e-08,
          substance=Chemical.Examples.Substances.PhosphoricAcid_aqueous)
          annotation (Placement(transformation(extent={{-98,-58},{-78,-38}})));
        Components.Substance H2PO4(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          solute_start=0.0005,
          substance=Chemical.Examples.Substances.DihydrogenPhosphate_aqueous)
          annotation (Placement(transformation(extent={{-44,-58},{-24,-38}})));
        Components.Substance HPO4(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          substance=Chemical.Examples.Substances.Phosphate[3],
          solute_start=0.0006)
          annotation (Placement(transformation(extent={{16,-58},{36,-38}})));
        Components.Substance PO4(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          substance=Chemical.Examples.Substances.Phosphate[4],
          solute_start=1e-08)
          annotation (Placement(transformation(extent={{72,-58},{92,-38}})));

        Chemical.Components.Reaction chemicalReaction(nP=2) "10^(-1.915 + 3)"
          annotation (Placement(transformation(extent={{-70,-58},{-50,-38}})));
        Chemical.Components.Reaction chemicalReaction1(nP=2) "10^(-6.66 + 3)"
          annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
        Chemical.Components.Reaction chemicalReaction2(nP=2) "10^(-11.78 + 3)"
          annotation (Placement(transformation(extent={{44,-58},{64,-38}})));

      equation
        connect(H3PO4.q_out, chemicalReaction.substrates[1]) annotation (Line(
            points={{-88,-48},{-70,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction.products[1], H2PO4.q_out) annotation (Line(
            points={{-50,-48.5},{-42,-48.5},{-42,-48},{-34,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2PO4.q_out, chemicalReaction1.substrates[1]) annotation (Line(
            points={{-34,-48},{-14,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[1], HPO4.q_out) annotation (Line(
            points={{6,-48.5},{16,-48.5},{16,-48},{26,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(HPO4.q_out, chemicalReaction2.substrates[1]) annotation (Line(
            points={{26,-48},{44,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction2.products[1], PO4.q_out) annotation (Line(
            points={{64,-48.5},{74,-48.5},{74,-48},{82,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction.products[2], H.q_out) annotation (Line(
            points={{-50,-47.5},{-44,-47.5},{-44,-32},{36,-32},{36,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[2], H.q_out) annotation (Line(
            points={{6,-47.5},{14,-47.5},{14,-32},{36,-32},{36,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction2.products[2], H.q_out) annotation (Line(
            points={{64,-47.5},{72,-47.5},{72,-32},{36,-32},{36,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>Henderson-Hasselbalch equation in ideal buffered solution, where pH remains constant.</p>
<p>The partial pressure of CO2 in gas are input parameter. Outputs are an amount of free disolved CO2 in liquid and an amount of HCO3-.</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.05));
      end Phosphate;

      model AlbuminTitration "Figge-Fencl model (22. Dec. 2007)"
        extends Modelica.Icons.Example;

        Physiolibrary.Chemical.Components.Substance H(
          q_out(conc(nominal=10^(-7.4 + 3))),
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          solute_start=10^(-7.4 + 3),
          isDependent=true) "hydrogen ions activity" annotation (Placement(
              transformation(extent={{-10,-10},{10,10}}, origin={14,22})));

        Physiolibrary.SteadyStates.Components.MolarConservationLaw molarConservationLaw[n](
          each n=2,
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          each Total=0.00066)
          annotation (Placement(transformation(extent={{44,-6},{64,14}})));
        Physiolibrary.SteadyStates.Components.ElementaryChargeConservationLaw electroneutrality(
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          NumberOfParticles=n,
          Charges=ones(n),
          useTotalInput=true,
          Total=6425.92363734) "strong ion difference of solution"
          annotation (Placement(transformation(extent={{46,-94},{66,-74}})));
        Modelica.Blocks.Math.Gain toColoumn(k(unit="C/s")=-Modelica.Constants.F,y(unit="C"))
          "from elementary charge to to electric charge, which is needed in system"
                                              annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={80,-74})));
        Modelica.Blocks.Sources.Clock SID(offset=-0.0832)
          "strong ions difference with respect to albumin charge shift"
          annotation (Placement(transformation(extent={{54,76},{74,96}})));

        parameter Integer n=218 "Number of weak acid group in albumin molecule";
        parameter Real pKAs[n]=cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
          "acid dissociation constants";

        Physiolibrary.Chemical.Components.Substance A[n](
          each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          each isDependent=true,
          each solute_start=0.00033) "deprotonated acid groups"
          annotation (Placement(transformation(extent={{4,-16},{24,4}})));
        Physiolibrary.Chemical.Components.ChemicalReaction react[n](each nP=2,
            K=fill(10.0, n) .^ (-pKAs .+ 3))
          annotation (Placement(transformation(extent={{-44,-2},{-24,18}})));

        Physiolibrary.Chemical.Components.Substance HA[n](each Simulation=
              Physiolibrary.Types.SimulationType.SteadyState, each solute_start=
             0.00033) "protonated acid groups"
          annotation (Placement(transformation(extent={{-76,-2},{-56,18}})));

      equation
        connect(react.products[1], A.q_out) annotation (Line(
            points={{-24,7.5},{-12,7.5},{-12,-6},{14,-6}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        for i in 1:n loop
          connect(react[i].products[2], H.q_out) annotation (Line(
              points={{-24,8.5},{-14,8.5},{-14,22},{14,22}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
        end for;
        connect(HA.q_out, react.substrates[1]) annotation (Line(
            points={{-66,8},{-44,8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(A.solute, molarConservationLaw.fragment[1]) annotation (Line(
            points={{20,-16},{20,-20},{36,-20},{36,-1},{44,-1}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(HA.solute, molarConservationLaw.fragment[2]) annotation (Line(
            points={{-60,-2},{-60,-8},{-78,-8},{-78,36},{36,36},{36,0},{44,0},{
                44,1}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(SID.y,toColoumn. u) annotation (Line(
            points={{75,86},{100,86},{100,-74},{92,-74}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(HA.solute, electroneutrality.fragment) annotation (Line(
            points={{-60,-2},{-60,-88},{46,-88}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(toColoumn.y, electroneutrality.total) annotation (Line(
            points={{69,-74},{56,-74},{56,-76}},
            color={0,0,127},
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
            StopTime=0.0235,
            __Dymola_fixedstepsize=5e-005,
            __Dymola_Algorithm="Euler"));
      end AlbuminTitration;

      class Develop
        extends Modelica.Icons.UnderConstruction;
        model PlasmaAcidBase

          Physiolibrary.Chemical.Components.Substance H3O(
            q_out(conc(nominal=10^(-7.4 + 3))),
            Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            solute_start=10^(-7.4 + 3),
            isDependent=true) "hydrogen ions activity" annotation (Placement(
                transformation(extent={{-10,-10},{10,10}}, origin={38,40})));

          Physiolibrary.SteadyStates.Components.MolarConservationLaw tAlb[n](
            each n=2,
            each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            each Total=0.00066)
            annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
          Physiolibrary.SteadyStates.Components.ElementaryChargeConservationLaw
            electroneutrality(
            Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            useTotalInput=true,
            Charges=cat(1,
                        {-1,-1,-2,-1},
                        ones(n)),
            NumberOfParticles=m + n,
            Total=6425.92363734) "strong ion difference of solution"
            annotation (Placement(transformation(extent={{46,-94},{66,-74}})));
          Modelica.Blocks.Math.Gain toColoumn(k(unit="C/s")=-Modelica.Constants.F,y(unit="C"))
            "from elementary charge to to electric charge, which is needed in system"
                                                annotation (Placement(transformation(
                extent={{-8,-8},{8,8}},
                rotation=180,
                origin={78,-70})));
          Modelica.Blocks.Sources.Clock SID_less_Cl(offset=-0.0832)
            "strong ions difference without chloride with respect to albumin charge shift"
            annotation (Placement(transformation(extent={{68,-42},{88,-22}})));

          constant Integer m=4
            "number of particle types in electroneutrality equation";

          parameter Boolean isDependent[3] = {false,false,false};

          parameter Physiolibrary.Types.AmountOfSubstance totalPO4=0.00115
            "Total phosphate concentration";
          parameter Physiolibrary.Types.AmountOfSubstance totalAlb=0.00066
            "Total albumin concentration";

          parameter Integer n=218
            "Number of weak acid group in albumin molecule";
          parameter Real pKAs[n]=cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
            "acid dissociation constants";

          Physiolibrary.Chemical.Components.Substance A[n](
            each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            each isDependent=true,
            each solute_start=0.00033) "deprotonated acid groups"
            annotation (Placement(transformation(extent={{-10,14},{10,34}})));
          Physiolibrary.Chemical.Components.ChemicalReaction react[n](each nP=2,
              K=fill(10.0, n) .^ (-pKAs .+ 3))
            annotation (Placement(transformation(extent={{-44,16},{-24,36}})));

          Physiolibrary.Chemical.Components.Substance HA[n](each Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, each
              solute_start=0.00033) "protonated acid groups"
            annotation (Placement(transformation(extent={{-76,16},{-56,36}})));

          Physiolibrary.Chemical.Components.Substance CO2_liquid(Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, isDependent=
                isDependent[1])
            annotation (Placement(transformation(extent={{-76,64},{-56,84}})));
          Physiolibrary.Chemical.Components.Substance HCO3(Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, isDependent=
                isDependent[2])
            annotation (Placement(transformation(extent={{42,70},{62,90}})));
          Physiolibrary.Chemical.Interfaces.ChemicalPort_a substances[3]
            "{free dissolved CO2, bicarbonate, chloride}"
            annotation (Placement(transformation(extent={{-10,70},{10,90}})));
          Physiolibrary.Chemical.Components.Substance H2PO4(Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, solute_start=
                0.0005) annotation (Placement(transformation(extent={{-62,-54},
                    {-42,-34}})));
          Physiolibrary.Chemical.Components.ChemicalReaction phosphateAcidification(nP=2, K=
                10^(-6.66 + 3)) annotation (Placement(transformation(extent={{-32,
                    -54},{-12,-34}})));
          Physiolibrary.Chemical.Components.Substance HPO4(
            Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            isDependent=true,
            solute_start=0.0006)
            annotation (Placement(transformation(extent={{-2,-54},{18,-34}})));
          Physiolibrary.SteadyStates.Components.MolarConservationLaw tP04(
            each Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            each n=2,
            each Total=totalPO4)
            annotation (Placement(transformation(extent={{-28,-80},{-8,-60}})));
          Physiolibrary.Chemical.Components.Substance Cl(Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, isDependent=
                isDependent[3]) "chloride anion"
            annotation (Placement(transformation(extent={{76,42},{96,62}})));
        equation
          connect(react.products[1], A.q_out) annotation (Line(
              points={{-24,25.5},{-12,25.5},{-12,24},{0,24}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          for i in 1:n loop
            connect(react[i].products[2], H3O.q_out) annotation (Line(
                points={{-24,26.5},{-14,26.5},{-14,40},{38,40}},
                color={107,45,134},
                thickness=1,
                smooth=Smooth.None));
          end for;
          connect(HA.q_out, react.substrates[1]) annotation (Line(
              points={{-66,26},{-44,26}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(A.solute, tAlb.fragment[1]) annotation (Line(
              points={{6,14},{6,10},{-54,10},{-54,-5},{-40,-5}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(HA.solute, tAlb.fragment[2]) annotation (Line(
              points={{-60,16},{-60,-4},{-40,-4},{-40,-3}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(toColoumn.y, electroneutrality.total) annotation (Line(
              points={{69.2,-70},{56,-70},{56,-76}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(substances[1], CO2_liquid.q_out) annotation (Line(
              points={{0,73.3333},{0,74},{-66,74}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(HCO3.q_out, substances[2]) annotation (Line(
              points={{52,80},{0,80}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(H2PO4.q_out, phosphateAcidification.substrates[1]) annotation (Line(
              points={{-52,-44},{-32,-44}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(phosphateAcidification.products[1], HPO4.q_out) annotation (Line(
              points={{-12,-44.5},{-2,-44.5},{-2,-44},{8,-44}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(H2PO4.solute, tP04.fragment[1]) annotation (Line(
              points={{-46,-54},{-46,-75},{-28,-75}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(HPO4.solute, tP04.fragment[2]) annotation (Line(
              points={{14,-54},{14,-60},{-40,-60},{-40,-73},{-28,-73}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(phosphateAcidification.products[2], H3O.q_out) annotation (Line(
              points={{-12,-43.5},{-4,-43.5},{-4,-28},{20,-28},{20,40},{38,40}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(HCO3.solute, electroneutrality.fragment[1]) annotation (Line(
              points={{58,70},{58,-62},{32,-62},{32,-88},{46,-88}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(H2PO4.solute, electroneutrality.fragment[2]) annotation (Line(
              points={{-46,-54},{-46,-86},{46,-86},{46,-88}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(HPO4.solute, electroneutrality.fragment[3]) annotation (Line(
              points={{14,-54},{14,-88},{46,-88}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(Cl.q_out, substances[3]) annotation (Line(
              points={{86,52},{0,52},{0,86.6667}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(HA.solute, electroneutrality.fragment[(m+1):(n+m)]) annotation (Line(
              points={{-60,16},{-60,-88},{46,-88}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(SID_less_Cl.y, toColoumn.u) annotation (Line(
              points={{89,-32},{92,-32},{92,-70},{87.6,-70}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(Cl.solute, electroneutrality.fragment[4]) annotation (Line(
              points={{92,42},{92,10},{54,10},{54,-64},{34,-64},{34,-88},{46,
                  -88}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation ( Documentation(revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",         info="<html>
<pre><b>plotExpression(apply(-log10(AlbuminTitration.H3O.solute)),&nbsp;false,&nbsp;&QUOT;pH&QUOT;,&nbsp;1);</b></pre>
<p>The titration slope der(pH)/der(SID)=185 1/(mol/L) at pH=7.4 and tAlb=0.66 mmol/l.</p>
<p><br>Data and model is described in</p>
<p><font style=\"color: #222222; \">Jame Figge: Role of non-volatile weak acids (albumin, phosphate and citrate). In: Stewart&apos;s Textbook of Acid-Base, 2nd Edition, John A. Kellum, Paul WG Elbers editors, &nbsp;AcidBase org, 2009, pp. 216-232.</font></p>
</html>"),  experiment(
              StopTime=0.0235,
              __Dymola_fixedstepsize=5e-005,
              __Dymola_Algorithm="Euler"));
        end PlasmaAcidBase;

        model ErythrocyteAcidBase
          parameter Boolean isDependent[4] = {false,false,false,false};

          Physiolibrary.Chemical.Components.Substance H3O(
            q_out(conc(nominal=10^(-7.4 + 3))),
            Simulation=Physiolibrary.Types.SimulationType.SteadyState,
            solute_start=10^(-7.4 + 3),
            isDependent=isDependent[4]) "hydrogen ions activity" annotation (
              Placement(transformation(extent={{-10,-10},{10,10}}, origin={-12,
                    36})));
          Physiolibrary.Chemical.Components.ChemicalReaction HendersonHasselbalch(
            nP=2,
            dH=15.13,
            K=10^(-6.103 + 3),
            nS=1)
            annotation (Placement(transformation(extent={{-60,46},{-40,66}})));
          Physiolibrary.Chemical.Components.Substance CO2_liquid(Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, isDependent=
                isDependent[1])
            annotation (Placement(transformation(extent={{-90,46},{-70,66}})));
          Physiolibrary.Chemical.Components.Substance HCO3(Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, isDependent=
                isDependent[2])
            annotation (Placement(transformation(extent={{-22,70},{-2,90}})));
          Physiolibrary.Chemical.Interfaces.ChemicalPort_a substances[3]
            "{free dissolved CO2, bicarbonate, chloride}"
            annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
          Physiolibrary.Chemical.Components.Substance Cl(Simulation=
                Physiolibrary.Types.SimulationType.SteadyState, isDependent=
                isDependent[3]) "chloride anion"
            annotation (Placement(transformation(extent={{76,82},{96,102}})));
        equation
          connect(HendersonHasselbalch.products[1],HCO3. q_out) annotation (Line(
              points={{-40,55.5},{-30,55.5},{-30,80},{-12,80}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(H3O.q_out,HendersonHasselbalch. products[2]) annotation (Line(
              points={{-12,36},{-30,36},{-30,56.5},{-40,56.5}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(CO2_liquid.q_out,HendersonHasselbalch. substrates[1]) annotation (
             Line(
              points={{-80,56},{-60,56}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(substances[1], CO2_liquid.q_out) annotation (Line(
              points={{-80,73.3333},{-80,56}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(HCO3.q_out, substances[2]) annotation (Line(
              points={{-12,80},{-80,80}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(Cl.q_out, substances[3]) annotation (Line(
              points={{86,92},{-80,92},{-80,86.6667}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
        end ErythrocyteAcidBase;

        model BloodAcidBase
          PlasmaAcidBase plasmaAcidBase
            annotation (Placement(transformation(extent={{-34,-6},{-14,14}})));
          Develop.ErythrocyteAcidBase erythrocyteAcidBase
            annotation (Placement(transformation(extent={{56,-6},{76,14}})));
          Physiolibrary.Chemical.Components.Membrane membrane(NumberOfParticles=
               3)
            annotation (Placement(transformation(extent={{14,2},{34,22}})));
        equation
          connect(plasmaAcidBase.substances, membrane.particlesInside)
            annotation (Line(
              points={{-24,12},{14,12}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(membrane.particlesOutside, erythrocyteAcidBase.substances)
            annotation (Line(
              points={{34,12},{58,12}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
        end BloodAcidBase;
      end Develop;

      model CarbonDioxideInBlood
        import Chemical;
          extends Modelica.Icons.Example;
        Components.Substance HCO3(
          substance=Chemical.Examples.Substances.CarbonicAcid[2],
          SolutionAmount=52.3,
          solute_start=0.024)
          annotation (Placement(transformation(extent={{24,22},{44,42}})));
        Sources.AirSubstance CO2_gas(
          substance=Chemical.Examples.Substances.CarbonDioxide[1],
          PartialPressure=5999.507433675,
          TotalPressure=101325.0144354) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-80,82})));
        Components.GasSolubility gasSolubility(kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")=
               0.00062064026806947)
          annotation (Placement(transformation(extent={{-90,46},{-70,66}})));
                                              /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
        Components.Substance CO2_liquid(
          substance=Chemical.Examples.Substances.CarbonDioxide[2],
          SolutionAmount=52.3,
          solute_start=0.00123)
          annotation (Placement(transformation(extent={{-92,22},{-72,42}})));
        Components.Substance H2O(
          substance=Chemical.Examples.Substances.Water,
          SolutionAmount=52.3,
          solute_start=52)
          annotation (Placement(transformation(extent={{-64,6},{-44,26}})));
        Components.Membrane membrane(
          Charges={0,0,-1,-1},
          NumberOfParticles=4,
          Permeabilities=10*ones(4))
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-4,-26})));
        Components.Substance HCO3_E(
          substance=Chemical.Examples.Substances.CarbonicAcid[2],
          SolutionAmount=39.7,
          solute_start=0.0116)
          annotation (Placement(transformation(extent={{26,-70},{46,-50}})));
        Chemical.Components.Reaction HendersonHasselbalch1(nP=2, nS=2)
          "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-16,-70},{4,-50}})));
        Components.Substance CO2_liquid_E(
          substance=Chemical.Examples.Substances.CarbonDioxide[2],
          SolutionAmount=39.7,
          solute_start=0.00093)
          annotation (Placement(transformation(extent={{-90,-70},{-70,-50}})));
        Components.Substance H2O_E(
          substance=Chemical.Examples.Substances.Water,
          SolutionAmount=39.7,
          solute_start=39.5)
          annotation (Placement(transformation(extent={{-64,-88},{-44,-68}})));
        Components.Substance Cl_E(
          substance=Chemical.Examples.Substances.Chloride,
          SolutionAmount=39.7,
          solute_start=0.0499)
          annotation (Placement(transformation(extent={{78,-94},{98,-74}})));
        Components.Substance Cl_P(
          substance=Chemical.Examples.Substances.Chloride,
          SolutionAmount=52.3,
          solute_start=0.103)
          annotation (Placement(transformation(extent={{78,-2},{98,18}})));

        Real pH_e;
      //  Real pH_p;
        Chemical.Sources.AmbientMolality H_E(substance=Chemical.Examples.Substances.Proton,
            Concentration(displayUnit="1") = 10^(-7.2)) annotation (Placement(
              transformation(extent={{10,-10},{-10,10}}, origin={60,-82})));
      equation
      //  pH_p = -log10(H.q_out.conc);
        pH_e = -log10(H_E.q_out.conc);
        connect(gasSolubility.q_in, CO2_liquid.q_out) annotation (Line(
            points={{-80,48},{-80,32},{-82,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_gas.q_out, gasSolubility.q_out) annotation (Line(
            points={{-80,72},{-80,66}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch1.products[1], HCO3_E.q_out) annotation (Line(
            points={{4,-60.5},{16,-60.5},{16,-60},{36,-60}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid_E.q_out, HendersonHasselbalch1.substrates[1]) annotation (
            Line(
            points={{-80,-60},{-28,-60},{-28,-60.5},{-16,-60.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O_E.q_out, HendersonHasselbalch1.substrates[2]) annotation (Line(
            points={{-54,-78},{-22,-78},{-22,-59.5},{-16,-59.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid_E.q_out, membrane.particlesOutside[1]) annotation (Line(
            points={{-80,-60},{-80,-44},{-4.75,-44},{-4.75,-36}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O_E.q_out, membrane.particlesOutside[2]) annotation (Line(
            points={{-54,-78},{-54,-44},{-4,-44},{-4,-36},{-4.25,-36}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HCO3_E.q_out, membrane.particlesOutside[3]) annotation (Line(
            points={{36,-60},{36,-44},{-3.75,-44},{-3.75,-36}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.q_out, membrane.particlesInside[1]) annotation (Line(
            points={{-82,32},{-82,-10},{-4.75,-10},{-4.75,-16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.q_out, membrane.particlesInside[2]) annotation (Line(
            points={{-54,16},{-54,-10},{-4.25,-10},{-4.25,-16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HCO3.q_out, membrane.particlesInside[3]) annotation (Line(
            points={{34,32},{34,-10},{-3.75,-10},{-3.75,-16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(Cl_P.q_out, membrane.particlesInside[4]) annotation (Line(
            points={{88,8},{88,-10},{-3.25,-10},{-3.25,-16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(Cl_E.q_out, membrane.particlesOutside[4]) annotation (Line(
            points={{88,-84},{88,-44},{-3.25,-44},{-3.25,-36}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H_E.q_out, HendersonHasselbalch1.products[2]) annotation (Line(
            points={{50,-82},{20,-82},{20,-59.5},{4,-59.5}},
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

      model CarbonDioxideInBlood2
        import Chemical;
          extends Modelica.Icons.Example;
        Sources.AirSubstance CO2_gas(
          substance=Chemical.Examples.Substances.CarbonDioxide[1],
          PartialPressure=5332.8954966,
          TotalPressure=101325.0144354) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-80,82})));
        Components.GasSolubility gasSolubility(kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")=
               0.00062064026806947)
          annotation (Placement(transformation(extent={{-90,46},{-70,66}})));
                                              /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
        Components.Substance CO2_liquid(
          substance=Chemical.Examples.Substances.CarbonDioxide[2],
          SolutionAmount=52.3,
          solute_start=0.00123)
          annotation (Placement(transformation(extent={{-92,22},{-72,42}})));
        Components.Substance HCO3(
          substance=Chemical.Examples.Substances.CarbonicAcid[2],
          SolutionAmount=52.3,
          solute_start=0.024)
          annotation (Placement(transformation(extent={{32,22},{52,42}})));
        Chemical.Components.Reaction HendersonHasselbalch(nP=2, nS=2)
          "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-12,22},{8,42}})));
        Components.Substance H(
          substance=Chemical.Examples.Substances.Proton,
          q_out(conc(nominal=10^(-7.4))),
          solute_start=52.3*10^(-7.4),
          SolutionAmount=52.3) annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={66,8})));
        Components.Substance H2O(
          substance=Chemical.Examples.Substances.Water,
          SolutionAmount=52.3,
          solute_start=52)
          annotation (Placement(transformation(extent={{-56,6},{-36,26}})));

        Real pH;
      equation
        pH=-log10(H.q_out.conc);
        connect(gasSolubility.q_in, CO2_liquid.q_out) annotation (Line(
            points={{-80,48},{-80,32},{-82,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_gas.q_out, gasSolubility.q_out) annotation (Line(
            points={{-80,72},{-80,66}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch.products[1],HCO3. q_out) annotation (Line(
            points={{8,31.5},{14,31.5},{14,32},{42,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H.q_out,HendersonHasselbalch. products[2]) annotation (Line(
            points={{66,8},{20,8},{20,32.5},{8,32.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.q_out,HendersonHasselbalch. substrates[1]) annotation (
           Line(
            points={{-82,32},{-62,32},{-62,31.5},{-12,31.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.q_out,HendersonHasselbalch. substrates[2]) annotation (Line(
            points={{-46,16},{-26,16},{-26,32.5},{-12,32.5}},
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
</html>"),experiment(__Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics),
          __Dymola_experimentSetupOutput);
      end CarbonDioxideInBlood2;

      model CarbonDioxideInBlood3
        import Chemical;
          extends Modelica.Icons.Example;
        Components.Substance HCO3(
          substance=Chemical.Examples.Substances.CarbonicAcid[2],
          SolutionAmount=52.3,
          solute_start=0.024)
          annotation (Placement(transformation(extent={{24,22},{44,42}})));
        Chemical.Components.Reaction HendersonHasselbalch(nP=2, nS=2)
          "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-20,22},{0,42}})));
        Sources.AirSubstance CO2_gas(
          substance=Chemical.Examples.Substances.CarbonDioxide[1],
          PartialPressure=5332.8954966,
          TotalPressure=101325.0144354) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-80,82})));
        Chemical.Sources.AmbientMolality H(substance=Chemical.Examples.Substances.Proton,
            Concentration(displayUnit="1") = 10^(-7.4)) annotation (Placement(
              transformation(extent={{-10,-10},{10,10}}, origin={58,8})));
        Components.GasSolubility gasSolubility(kH_T0(displayUnit=
                "(mol/kg H2O)/bar at 25degC,101325Pa") = 0.00062064026806947)
          annotation (Placement(transformation(extent={{-90,46},{-70,66}})));
                                              /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
        Components.Substance CO2_liquid(
          substance=Chemical.Examples.Substances.CarbonDioxide[2],
          SolutionAmount=52.3,
          solute_start=0.00123)
          annotation (Placement(transformation(extent={{-92,22},{-72,42}})));
        Components.Substance H2O(
          substance=Chemical.Examples.Substances.Water,
          SolutionAmount=52.3,
          solute_start=52)
          annotation (Placement(transformation(extent={{-64,6},{-44,26}})));
        Components.Membrane membrane(
          Charges={0,0,-1,-1},
          NumberOfParticles=4,
          Permeabilities=10*ones(4))
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-4,-26})));
        Components.Substance HCO3_E(
          substance=Chemical.Examples.Substances.CarbonicAcid[2],
          SolutionAmount=39.7,
          solute_start=0.0116)
          annotation (Placement(transformation(extent={{26,-70},{46,-50}})));
        Chemical.Components.Reaction HendersonHasselbalch1(nP=2, nS=2)
          "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-16,-70},{4,-50}})));
        Components.Substance H_E(
          substance=Chemical.Examples.Substances.Proton,
          solute_start=39.7*10^(-7.2),
          SolutionAmount=39.7) annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={60,-84})));
        Components.Substance CO2_liquid_E(
          substance=Chemical.Examples.Substances.CarbonDioxide[2],
          SolutionAmount=39.7,
          solute_start=0.00093)
          annotation (Placement(transformation(extent={{-90,-70},{-70,-50}})));
        Components.Substance H2O_E(
          substance=Chemical.Examples.Substances.Water,
          SolutionAmount=39.7,
          solute_start=39.5)
          annotation (Placement(transformation(extent={{-64,-88},{-44,-68}})));
        Components.Substance Cl_E(
          substance=Chemical.Examples.Substances.Chloride,
          SolutionAmount=39.7,
          solute_start=0.0499)
          annotation (Placement(transformation(extent={{78,-94},{98,-74}})));
        Components.Substance Cl_P(
          substance=Chemical.Examples.Substances.Chloride,
          SolutionAmount=52.3,
          solute_start=0.103)
          annotation (Placement(transformation(extent={{78,-2},{98,18}})));

        Real pH_e,pH_p;
      equation
        pH_p = -log10(H.q_out.conc);
        pH_e = -log10(H_E.q_out.conc);
        connect(HendersonHasselbalch.products[1], HCO3.q_out) annotation (Line(
            points={{0,31.5},{6,31.5},{6,32},{34,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H.q_out, HendersonHasselbalch.products[2]) annotation (Line(
            points={{68,8},{12,8},{12,32.5},{0,32.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.q_out, HendersonHasselbalch.substrates[1]) annotation (
           Line(
            points={{-82,32},{-70,32},{-70,31.5},{-20,31.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(gasSolubility.q_in, CO2_liquid.q_out) annotation (Line(
            points={{-80,48},{-80,32},{-82,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_gas.q_out, gasSolubility.q_out) annotation (Line(
            points={{-80,72},{-80,66}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.q_out, HendersonHasselbalch.substrates[2]) annotation (Line(
            points={{-54,16},{-34,16},{-34,32.5},{-20,32.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch1.products[1], HCO3_E.q_out) annotation (Line(
            points={{4,-60.5},{16,-60.5},{16,-60},{36,-60}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H_E.q_out, HendersonHasselbalch1.products[2]) annotation (Line(
            points={{60,-84},{16,-84},{16,-59.5},{4,-59.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid_E.q_out, HendersonHasselbalch1.substrates[1]) annotation (
            Line(
            points={{-80,-60},{-28,-60},{-28,-60.5},{-16,-60.5}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O_E.q_out, HendersonHasselbalch1.substrates[2]) annotation (Line(
            points={{-54,-78},{-22,-78},{-22,-59.5},{-16,-59.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid_E.q_out, membrane.particlesOutside[1]) annotation (Line(
            points={{-80,-60},{-80,-44},{-4.75,-44},{-4.75,-36}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O_E.q_out, membrane.particlesOutside[2]) annotation (Line(
            points={{-54,-78},{-54,-44},{-4,-44},{-4,-36},{-4.25,-36}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HCO3_E.q_out, membrane.particlesOutside[3]) annotation (Line(
            points={{36,-60},{36,-44},{-3.75,-44},{-3.75,-36}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.q_out, membrane.particlesInside[1]) annotation (Line(
            points={{-82,32},{-82,-10},{-4.75,-10},{-4.75,-16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.q_out, membrane.particlesInside[2]) annotation (Line(
            points={{-54,16},{-54,-10},{-4.25,-10},{-4.25,-16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HCO3.q_out, membrane.particlesInside[3]) annotation (Line(
            points={{34,32},{34,-10},{-3.75,-10},{-3.75,-16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(Cl_P.q_out, membrane.particlesInside[4]) annotation (Line(
            points={{88,8},{88,-10},{-3.25,-10},{-3.25,-16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(Cl_E.q_out, membrane.particlesOutside[4]) annotation (Line(
            points={{88,-84},{88,-44},{-3.25,-44},{-3.25,-36}},
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
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}),            graphics),
          __Dymola_experimentSetupOutput);
      end CarbonDioxideInBlood3;

      model Phosphate2
        import Chemical;
          extends Modelica.Icons.Example;

        parameter Physiolibrary.Types.Concentration totalPO4=0.00115
          "Total phosphate concentration";

        Modelica.Blocks.Sources.Clock SID(offset=-1.9*totalPO4)
          "strong ions difference with respect to albumin charge shift"
          annotation (Placement(transformation(extent={{44,74},{64,94}})));

        Components.Substance H(
          solute_start=55.6*10^(-7.4),
          q_out(conc(nominal=10^(-7.4))),
          substance=Chemical.Examples.Substances.Proton,
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isDependent=true,
          NominalSolute(displayUnit="mol") = 1e-06) "hydrogen ions activity"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                origin={36,-12})));

        Components.Substance H3PO4(
          substance=Chemical.Examples.Substances.Phosphate[1],
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          solute_start=4.8e-10,
          NominalSolute=1e-10)
          annotation (Placement(transformation(extent={{-98,-58},{-78,-38}})));
        Components.Substance H2PO4(
          substance=Chemical.Examples.Substances.Phosphate[2],
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          solute_start=0.0005)
          annotation (Placement(transformation(extent={{-44,-58},{-24,-38}})));
        Components.Substance HPO4(
          substance=Chemical.Examples.Substances.Phosphate[3],
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isDependent=true,
          solute_start=0.0006)
          annotation (Placement(transformation(extent={{16,-58},{36,-38}})));
        Components.Substance PO4(
          substance=Chemical.Examples.Substances.Phosphate[4],
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          solute_start=2.3e-07,
          NominalSolute=1e-07)
          annotation (Placement(transformation(extent={{72,-58},{92,-38}})));

        Chemical.Components.Reaction chemicalReaction(nP=2) "10^(-1.915 + 3)"
          annotation (Placement(transformation(extent={{-70,-58},{-50,-38}})));
        Chemical.Components.Reaction chemicalReaction1(nP=2) "10^(-6.66 + 3)"
          annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
        Chemical.Components.Reaction chemicalReaction2(nP=2) "10^(-11.78 + 3)"
          annotation (Placement(transformation(extent={{44,-58},{64,-38}})));
        Physiolibrary.SteadyStates.Components.MolarConservationLaw tP04(
          each n=4,
          each Total=totalPO4*1,
          Simulation=Physiolibrary.Types.SimulationType.SteadyState)
          annotation (Placement(transformation(extent={{-28,-90},{-8,-70}})));

        Modelica.Blocks.Math.Gain toColoumn(k(unit="C/s")=Modelica.Constants.F,y(unit="C"))
          "from elementary charge to Coloumn" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={84,-86})));

        Physiolibrary.SteadyStates.Components.ElementaryChargeConservationLaw electroneutrality(
          Total(displayUnit="meq") = 3502.41783837,
          useTotalInput=true,
          NumberOfParticles=3,
          Charges={-1,-2,-3},
          Simulation=Physiolibrary.Types.SimulationType.SteadyState)
          annotation (Placement(transformation(extent={{44,-104},{64,-84}})));
      equation
        connect(H3PO4.q_out, chemicalReaction.substrates[1]) annotation (Line(
            points={{-88,-48},{-70,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction.products[1], H2PO4.q_out) annotation (Line(
            points={{-50,-48.5},{-42,-48.5},{-42,-48},{-34,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2PO4.q_out, chemicalReaction1.substrates[1]) annotation (Line(
            points={{-34,-48},{-14,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[1], HPO4.q_out) annotation (Line(
            points={{6,-48.5},{16,-48.5},{16,-48},{26,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(HPO4.q_out, chemicalReaction2.substrates[1]) annotation (Line(
            points={{26,-48},{44,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction2.products[1], PO4.q_out) annotation (Line(
            points={{64,-48.5},{74,-48.5},{74,-48},{82,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction.products[2], H.q_out) annotation (Line(
            points={{-50,-47.5},{-44,-47.5},{-44,-32},{36,-32},{36,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[2], H.q_out) annotation (Line(
            points={{6,-47.5},{14,-47.5},{14,-32},{36,-32},{36,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction2.products[2], H.q_out) annotation (Line(
            points={{64,-47.5},{72,-47.5},{72,-32},{36,-32},{36,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H3PO4.solute, tP04.fragment[1]) annotation (Line(
            points={{-82,-58},{-82,-86},{-28,-86},{-28,-85.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(H2PO4.solute, tP04.fragment[2]) annotation (Line(
            points={{-28,-58},{-28,-62},{-64,-62},{-64,-84.5},{-28,-84.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(HPO4.solute, tP04.fragment[3]) annotation (Line(
            points={{32,-58},{32,-64},{-50,-64},{-50,-83.5},{-28,-83.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(PO4.solute, tP04.fragment[4]) annotation (Line(
            points={{88,-58},{88,-68},{-40,-68},{-40,-82.5},{-28,-82.5}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(SID.y, toColoumn.u) annotation (Line(
            points={{65,84},{100,84},{100,-86},{96,-86}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(H2PO4.solute,electroneutrality. fragment[1]) annotation (Line(
            points={{-28,-58},{-28,-62},{20,-62},{20,-99.3333},{44,-99.3333}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(HPO4.solute,electroneutrality. fragment[2]) annotation (Line(
            points={{32,-58},{32,-98},{44,-98}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(PO4.solute,electroneutrality. fragment[3]) annotation (Line(
            points={{88,-58},{88,-62},{24,-62},{24,-96.6667},{44,-96.6667}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(toColoumn.y, electroneutrality.total) annotation (Line(
            points={{73,-86},{68,-86},{68,-76},{54,-76},{54,-86}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>Henderson-Hasselbalch equation in ideal buffered solution, where pH remains constant.</p>
<p>The partial pressure of CO2 in gas are input parameter. Outputs are an amount of free disolved CO2 in liquid and an amount of HCO3-.</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(
            StopTime=0.0013,
            Tolerance=1e-006,
            __Dymola_Algorithm="Euler"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics),
          __Dymola_experimentSetupOutput);
      end Phosphate2;

      model Phosphate3
        import Chemical;
          extends Modelica.Icons.Example;

        parameter Physiolibrary.Types.Concentration totalPO4=0.00115
          "Total phosphate concentration";

        Modelica.Blocks.Sources.Clock SID(offset=-1.99*totalPO4)
          "strong ions difference with respect to albumin charge shift"
          annotation (Placement(transformation(extent={{44,74},{64,94}})));

        Components.Substance H(
          solute_start=55.6*10^(-7.4),
          q_out(conc(nominal=10^(-7.4))),
          substance=Chemical.Examples.Substances.Proton,
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isDependent=true,
          NominalSolute(displayUnit="mol") = 1e-06) "hydrogen ions activity"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                origin={36,-12})));

        Components.Substance H2PO4(
          substance=Chemical.Examples.Substances.Phosphate[2],
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          isDependent=true,
          solute_start=0.0005)
          annotation (Placement(transformation(extent={{-48,-54},{-28,-34}})));
        Components.Substance HPO4(
          substance=Chemical.Examples.Substances.Phosphate[3],
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          solute_start=0.0006)
          annotation (Placement(transformation(extent={{16,-58},{36,-38}})));

        Chemical.Components.Reaction chemicalReaction1(nP=2) "10^(-6.66 + 3)"
          annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
        Physiolibrary.SteadyStates.Components.MolarConservationLaw tP04(
          each Total=totalPO4*1,
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          each n=2)
          annotation (Placement(transformation(extent={{-28,-90},{-8,-70}})));

        Modelica.Blocks.Math.Gain toColoumn(k(unit="C/s")=Modelica.Constants.F,y(unit="C"))
          "from elementary charge to Coloumn" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={84,-86})));

        Physiolibrary.SteadyStates.Components.ElementaryChargeConservationLaw electroneutrality(
          Total(displayUnit="meq") = 3502.41783837,
          useTotalInput=true,
          Simulation=Physiolibrary.Types.SimulationType.SteadyState,
          NumberOfParticles=2,
          Charges={-1,-2})
          annotation (Placement(transformation(extent={{44,-104},{64,-84}})));
      equation
        connect(H2PO4.q_out, chemicalReaction1.substrates[1]) annotation (Line(
            points={{-38,-44},{-24,-44},{-24,-48},{-14,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[1], HPO4.q_out) annotation (Line(
            points={{6,-48.5},{16,-48.5},{16,-48},{26,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[2], H.q_out) annotation (Line(
            points={{6,-47.5},{14,-47.5},{14,-32},{36,-32},{36,-12}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(SID.y, toColoumn.u) annotation (Line(
            points={{65,84},{100,84},{100,-86},{96,-86}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(H2PO4.solute,electroneutrality. fragment[1]) annotation (Line(
            points={{-32,-54},{-32,-62},{20,-62},{20,-99},{44,-99}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(HPO4.solute,electroneutrality. fragment[2]) annotation (Line(
            points={{32,-58},{32,-97},{44,-97}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(toColoumn.y, electroneutrality.total) annotation (Line(
            points={{73,-86},{68,-86},{68,-76},{54,-76},{54,-86}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(H2PO4.solute, tP04.fragment[1]) annotation (Line(
            points={{-32,-54},{-54,-54},{-54,-85},{-28,-85}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(HPO4.solute, tP04.fragment[2]) annotation (Line(
            points={{32,-58},{-40,-58},{-40,-83},{-28,-83}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>Henderson-Hasselbalch equation in ideal buffered solution, where pH remains constant.</p>
<p>The partial pressure of CO2 in gas are input parameter. Outputs are an amount of free disolved CO2 in liquid and an amount of HCO3-.</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.001),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),
                          graphics),
          __Dymola_experimentSetupOutput);
      end Phosphate3;
    end AcidBase;

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

      parameter Boolean Isothermal = true
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

      if Isothermal then
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
            graphics={Text(
              extent={{-88,-92},{92,-100}},
              lineColor={0,0,255},
              textString="%name",
              horizontalAlignment=TextAlignment.Left),
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
    end Solution;

    model Substance "Substance in solution"
      extends Icons.Substance;

      parameter Modelica.SIunits.AmountOfSubstance amountOfSubstance_start=1e-8
        "Initial amount of the substance in compartment";

      replaceable package substance = Interfaces.BaseSubstance    constrainedby
        Interfaces.BaseSubstance
        "Substance definition: Molar Weight, Enthalpy, Gibbs energy,..."
         annotation (choicesAllMatching = true);

       Modelica.Blocks.Interfaces.RealOutput amountOfSubstance(start=amountOfSubstance_start, stateSelect=StateSelect.avoid, final unit="mol")
        "Current amount of the substance" annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,-60}),  iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,-60})));

      Interfaces.SubstanceDefinitionPort port_a
        "Concentration and molar flow from/to compartment" annotation (Placement(
            transformation(extent={{-10,-10},{10,10}}),  iconTransformation(extent={{-10,-10},
                {10,10}})));

      Interfaces.SolutionPort solution annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}),
            iconTransformation(extent={{-70,-110},{-50,-90}})));

    protected
      Modelica.SIunits.MoleFraction x(stateSelect=StateSelect.avoid)
        "Mole fraction of the substance";
      Modelica.SIunits.ChemicalPotential u
        "Electrochemical potential of the substance";

      Real state(stateSelect=StateSelect.prefer)
        "To help the numerical calculation: Using of logarithmic substitution for state variable as amountOfSubstance = exp(state), where der(amountOfSubstance) = der(state)*exp(state)";

      Modelica.SIunits.ChemicalPotential u0=
        substance.u0(solution.T,solution.p)
        "Base chemical potential of the pure substance";

      Modelica.SIunits.ChemicalPotential uPure=
        substance.uPure(solution.T,solution.v,solution.p)
        "Electro-chemical potential of the pure substance";

      Real activityCoefficient(final unit="1")=
          substance.activityCoefficient(solution.I)
        "Activity coefficient of the substance";

      Modelica.SIunits.MolarEnthalpy molarEnthalpy(displayUnit="kJ/mol")=
        substance.molarEnthalpy(solution.v,solution.p)
        "Molar enthalpy of the substance";

      Modelica.SIunits.MolarEntropy molarEntropy=
      substance.molarEntropy(port_a.u,solution.T,solution.v,solution.p)
        "Molar entropy of the substance";

      Modelica.SIunits.MolarVolume molarVolume=
      substance.molarVolume(solution.p,solution.T,solution.v,solution.I)
        "Molar volume of the substance";

    initial equation
      amountOfSubstance=amountOfSubstance_start;
    equation
      //der(amountOfSubstance)=port_a.q; <- This is mathematically the same as two following lines. However, the differential solvers can handle it much better. :-)
      der(state)=port_a.q/amountOfSubstance;
      amountOfSubstance = exp(state);

      //mole fraction (an analogy of molar concentration or molality)
      //if you select the amount of solution per one kilogram of solvent then the values of amountOfSubstance will be the same as molality
      //if you select the amount of solution in one liter of solution then the values of amountOfSubstance will be the same as molarity
      x = amountOfSubstance/solution.n;

      //electro-chemical potential of the substance in the solution
      u =
        u0
        + Modelica.Constants.R*solution.T*log(activityCoefficient*x)
        + substance.z*Modelica.Constants.F*solution.v;

      //define the substance at the port
      port_a.u = u;

      port_a.molarWeight = substance.MolarWeight;
      port_a.z = substance.z;
      port_a.x = x;
      port_a.molarEnthalpy = molarEnthalpy;
      port_a.molarEntropy = molarEntropy;
      port_a.molarVolume = molarVolume;
      port_a.u0 = u0;
      port_a.uPure = uPure;
      port_a.activityCoefficient = activityCoefficient;
      port_a.temperature = solution.T;
      port_a.electricPotential = solution.v;
      port_a.amountOfSolution = solution.n;

      //changes of the solution
      solution.dH = molarEnthalpy*port_a.q;
      solution.dS = molarEntropy*port_a.q;
      solution.dG = port_a.u*port_a.q;
      solution.dn = port_a.q;
      solution.i = Modelica.Constants.F * substance.z * port_a.q;
      solution.dI = (1/2) * port_a.q * substance.z^2;
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

      parameter Integer nP=1 "Number of products types"
        annotation ( HideResult=true);

      parameter Modelica.SIunits.Time Tau = 1
        "Time constant of kinetics: Tau * reaction molar rate = electro-chemical potential gradient * amount of solution";

      parameter Modelica.SIunits.StoichiometricNumber s[nS]=ones(nS)
        "Stoichiometric reaction coefficient for substrates"
        annotation (HideResult=true);

      parameter Modelica.SIunits.StoichiometricNumber p[nP]=ones(nP)
        "Stoichiometric reaction coefficients for products"
        annotation (HideResult=true);

      Modelica.SIunits.MolarFlowRate rr(start=0) "Reaction molar flow rate";

      Interfaces.SubstanceUsePort products[nP] "Products"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      Interfaces.SubstanceUsePort substrates[nS] "Substrates"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

      //solution properties:
      Modelica.SIunits.AmountOfSubstance amountOfSolution
        "amount of all particles in the solution";

      //for debugging:
      Real DissociationConstant
        "Dissociation constant as ratio of mole fractions";

      Modelica.SIunits.MolarEnergy DrH
        "Standard Enthalpy Change of reaction (negative=exothermic)";

      Modelica.SIunits.Power lossHeat "Comsumed heat by the reaction";

    //  Modelica.SIunits.ElectricPotential StandardNernstPotential
    //    "Standard electric potential of half-cell rection";

    protected
      Modelica.SIunits.MolarEnergy One=1 "Unit check";
    equation
      //the main equation
      rr = (1/(Tau*One))*amountOfSolution*((p * products.u) - (s * substrates.u));

      //reaction molar rates
      rr*s = -substrates.q;
      rr*p = products.q;

      //properties of solution
      amountOfSolution = substrates[1].amountOfSolution;

      //for debugging olny:

        //evaluation of dissociation constant from the Gibbs energy of the reaction
        DissociationConstant = (product(substrates.activityCoefficient .^ s) / product(products.activityCoefficient .^ p)) *
        exp( -(1/(Modelica.Constants.R*substrates[1].temperature)) * (p * products.uPure - s * substrates.uPure));

      //  0 = p * products.u0 - s * substrates.u0 + (p*products.z - s*substrates.z) * Modelica.Constants.F * StandardNernstPotential;

        //molar enthalpy of reaction
        DrH = sum(p.*products.molarEnthalpy) - sum(s.*substrates.molarEnthalpy);

        //consumed heat by reaction
        lossHeat = DrH*rr; //DrH<0 => Exothermic => lossHeat>0, Endothermic otherwise

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

       parameter Modelica.SIunits.Time Tau = 1
        "Time constant of kinetics: Tau * reaction molar rate = electro-chemical potential gradient";

    protected
      Real One(final unit="J/mol2")=1 "Unit check";
    equation
      port_b.q = (1/(Tau*One)) * (port_b.u - port_a.u);

       annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
    end Diffusion;

    model GasSolubility "Henry's law of gas solubility in liquid."

      extends Icons.GasSolubility;

      parameter Modelica.SIunits.Time Tau=1
        "Time constant of kinetics: Tau * dissolution molar rate = electro-chemical potential gradient";

      parameter Boolean useWaterCorrection = true
        "Are free Gibbs energy of aqueous formation shifted by 10 kJ/mol?"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));

      Interfaces.SubstanceUsePort gas_port "Gaseous solution"
        annotation (Placement(transformation(extent={{-10,90},{10,110}})));

      Interfaces.SubstanceUsePort liquid_port "Dissolved in liquid solution"
        annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
            iconTransformation(extent={{-10,-110},{10,-90}})));

      //for debugging
      Real kH(final unit="(mol/mol)/(Pa/Pa)", displayUnit="(mol/kg H2O)/bar at 25degC")
        "Henry's law coefficient such as liquid-gas concentration ratio at 25degC";
      Modelica.SIunits.Temperature C(displayUnit="K")
        "Henry's law temperature dependence coefficient";

      Modelica.SIunits.Power lossHeat "Comsumed heat by the reaction";

    //  Modelica.SIunits.ElectricPotential StandardNernstPotential
    //    "Standard electric potential";

    protected
      Real One(final unit="J/mol2")=1 "Unit check";
    equation
      gas_port.q + liquid_port.q = 0;

      // the main equation
      liquid_port.q = (1/(Tau*One))*(liquid_port.u - gas_port.u - (if useWaterCorrection then Modelica.Constants.R*liquid_port.temperature*log(0.018) else 0));

      lossHeat = (liquid_port.molarEnthalpy - gas_port.molarEnthalpy)*liquid_port.q; //negative = heat are comsumed when change from liquid to gas

      //for debugging olny:

        // evaluation of kH and C from enthalpies and Gibbs energies
        C=-(liquid_port.molarEnthalpy - gas_port.molarEnthalpy)/Modelica.Constants.R;

        -Modelica.Constants.R*liquid_port.temperature*
         log(kH*(liquid_port.activityCoefficient/gas_port.activityCoefficient)) =
         (liquid_port.uPure - gas_port.uPure)
         - (if useWaterCorrection then Modelica.Constants.R*liquid_port.temperature*log(0.018) else 0);

     //   0 = liquid_port.u0 - gas_port.u0 + (liquid_port.z - gas_port.z) * Modelica.Constants.F * StandardNernstPotential
     //    - (if useWaterCorrection then Modelica.Constants.R*liquid_port.temperature*log(0.018) else 0);

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

      parameter Modelica.SIunits.Time Tau = 1
        "Time constant of kinetics: Tau * membrane molar flux = electro-chemical potential gradient * membrane area";

      parameter Modelica.SIunits.Area MembraneArea=1e-4
        "Surface of the membrane";

      Interfaces.SubstanceUsePort port_a
        "The substance of inner side of membrane"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      Interfaces.SubstanceUsePort port_b
        "The substances of outer side of membrane"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

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

    protected
      Real One(final unit="J.m2/mol2")=1 "Unit check";
    equation
      //the main equation
      port_a.q = MembraneArea * (1 / (Tau*One)) * (port_a.u - port_b.u);

      //for debuging only:

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

      annotation ( Documentation(info="<html>
<p><u><b><font style=\"color: #008000; \">Filtration throught semipermeable membrane.</font></b></u></p>
<p>The penetrating particles are driven by electric and chemical gradient to reach Donnan&apos;s equilibrium.</p>
<p>If zero-flow Donnan&apos;s equilibrium is reached. </p>
</html>", revisions="<html>
<p><i>2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Membrane;

    model Speciation
      "Quaternary form of macromolecule with independent subunits"
      extends Icons.Speciation;

      replaceable package Macromolecule = Interfaces.BaseSubstance    constrainedby
        Interfaces.BaseSubstance
        "Substance definition: Molar Weight, Enthalpy, Gibbs energy,..."
        annotation (choicesAllMatching = true);

      parameter Integer NumberOfSubunits=1
        "Number of independent subunits occuring in macromolecule";

      Interfaces.SubstanceDefinitionPort macromolecule
        "Macromolecule composed with subunits"
        annotation (Placement(transformation(extent={{90,-90},{110,-70}})));
      Interfaces.SubstanceUsePort subunits[NumberOfSubunits]
        "Subunits of macromolecule"
        annotation (Placement(transformation(extent={{-10,90},{10,110}})));

      Modelica.Blocks.Interfaces.RealInput amountOfSubunits[
        NumberOfSubunits](each final unit="mol")
        "Total amount of the subunits in all their selected forms"                   annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=180,
            origin={80,0})));

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

      Modelica.SIunits.ActivityCoefficient activityCoefficient=
          Macromolecule.activityCoefficient(solution.I,solution.p)
        "Activity coefficient of the substance";

      Modelica.SIunits.MolarVolume molarVolume=
          Macromolecule.molarVolume(solution.p,solution.T,solution.v,solution.I)
        "Molar volume of the substance";

    protected
      Real fractions[NumberOfSubunits]
        "Fractions of selected specific form of each subunit in macromolecule";

    public
      Interfaces.SolutionPort solution annotation (Placement(transformation(extent={{-72,-90},{-52,-70}}),
            iconTransformation(extent={{-72,-90},{-52,-70}})));
    equation
      //change of macromolecule = change of its subunits
      subunits.q = -macromolecule.q * ones(NumberOfSubunits);

      //TODO: check if all subunits have the same amountOfSubunit as amountOfMacromolecule

      //the amount of total macromolecule is the same as amount of each its selected subunit
      amountOfMacromolecule = amountOfSubunits[1];

      //chemical speciation
      macromolecule.x = (amountOfMacromolecule/solution.n)*product(fractions);

      fractions = if (amountOfMacromolecule < Modelica.Constants.eps) then zeros(NumberOfSubunits)
                  else subunits.x ./ (amountOfSubunits/solution.n);

      macromolecule.u
      = Macromolecule.u0(solution.T,solution.p)
      + Modelica.Constants.R*solution.T*log(macromolecule.activityCoefficient .* macromolecule.x)
      + macromolecule.z * Modelica.Constants.F * solution.v;

      //define the macromolecule at the port
      macromolecule.activityCoefficient = activityCoefficient;
      macromolecule.molarWeight = Macromolecule.MolarWeight;
      macromolecule.z = Macromolecule.z;
      macromolecule.molarEnthalpy = Macromolecule.molarEnthalpy(solution.v,solution.p);
      macromolecule.molarEntropy = Macromolecule.molarEntropy(macromolecule.u,solution.T,solution.v,solution.p);
      macromolecule.molarVolume = molarVolume;
      macromolecule.u0 = Macromolecule.u0(solution.T,solution.p);
      macromolecule.uPure = Macromolecule.uPure(solution.T,solution.v,solution.p);
      macromolecule.temperature = solution.T;
      macromolecule.electricPotential = solution.v;
      macromolecule.amountOfSolution = solution.n;

      //changes of the solution, where all subunits are also connected in the same solution
      solution.dH = macromolecule.molarEnthalpy * macromolecule.q
      + subunits.molarEnthalpy * subunits.q;
      solution.dS = macromolecule.molarEntropy * macromolecule.q
      + subunits.molarEntropy * subunits.q;
      solution.dG = macromolecule.u * macromolecule.q + subunits.u * subunits.q;
      solution.dn = macromolecule.q + sum(subunits.q);
      solution.i = Modelica.Constants.F * (macromolecule.z * macromolecule.q + subunits.z * subunits.q);
      solution.dI = (1/2) * macromolecule.q * macromolecule.z^2
      + (1/2) * subunits.q * (subunits.z .^ 2);
      solution.dV = molarVolume * macromolecule.q
      + subunits.molarVolume * subunits.q;

      annotation (defaultComponentName="macromoleculeSpecie_in_macromoleculeGroup",
        Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>This block identifies one specific chemical form of one macromolecule defined by forms of its subunits  (one chosen chemical species - called <i>specie</i>).</p>
<p>Only main connector called <b>species </b>is designed for inflow and outflow of macromolecule to/from <i>system</i>. The concentration in this connector is the concentration of its specific <i>specie.</i></p>
<p>Connectors <b>subunitSpecies[:] </b>represent specific forms of the macromolecule subunit types. If the subnunit type occures n-times in macromolecule, the inflow is n-time greater than the inflow of macromolecule.</p>
<p><br>Initial total concentrations of subunits must be set to be right distribution of total macromolecule concentration. So the ratios between subunit concentrations are the ratios of their occurence in macromolecule. In equilibrium are this proporties fullfiled.</p>
<p><br>For example: If the macromolecule has four identical independent subunits and each subunit can occur in two form F1 and F2, then the concentration of macromolecule <i>specie </i>composed only from four subunits in form F1 is <b>species.conc=</b>conc*fF1^4. </p>
<p>Where:</p>
<p>conc is totat concentration of macromolecule in <i>system</i> accumulated by <b>species.q</b>,</p>
<p>fF1 = F1/(F1+F2) is fraction of form F1 in subsystem of subunit,</p>
<p>4 is number of subunits (<b>numberOfSubunit</b>).</p>
<p><br>This block can be connected to chemical reactions such as it was the chosen species with subsystem behind. It is recommended to use this block only as an equilibrated subsystem.</p>
<h4><span style=\"color:#008000\">Heat of chemical system.</span></h4>
<p>Enthalpy of each subunit species can be presented as <b>subunitEnthalpies[:]</b>. Then the total enthalpy of the chemical system can be calculated by equation:</p>
<h4>systemEnthalpy = sum(subunitEnthalpies[i] * totalSubunitAmount[i]) / totalSubsystemAmount</h4>
<p>And the stored heat as enthalpy is <b>systemEnthalpy*totalSubsystemAmount.</b></p>
</html>"),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics={                                                        Text(
              extent={{-22,-106},{220,-140}},
              lineColor={0,0,255},
              textString="%name")}));
    end Speciation;

    model Stream "Flow of whole solution"
      extends Interfaces.OnePortParallel;
      extends Interfaces.ConditionalSolutionFlow;

    equation
      port_a.q = if (q>0) then q*port_a.x else q*port_b.x;

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

      Interfaces.SubstanceUsePort port_a "For measure only"
        annotation (Placement(transformation(extent={{-10,-12},{10,8}}),
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

      moleFraction = port_a.x;

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

    parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionPer1kgOfSolvent = 55.508
        "Amount of all particles in the solution per one kilogram of solvent";

      Interfaces.SubstanceUsePort port_a "For measure only"
        annotation (Placement(transformation(extent={{-10,-12},{10,8}}),
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

      port_a.x=molality*KG / AmountOfSolutionPer1kgOfSolvent;

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

    parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionInOneLiter = 55.508
        "Amount of all particles in one liter of the solution";

      Interfaces.SubstanceUsePort port_a "For measure only"
        annotation (Placement(transformation(extent={{-10,-12},{10,8}}),
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

      port_a.x=molarConcentration*L / AmountOfSolutionInOneLiter;

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

    parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionInOneKilogram = 55.508
        "Amount of all particles in one kilogram of the solution";

      Interfaces.SubstanceUsePort port_a "For measure only"
        annotation (Placement(transformation(extent={{-10,-12},{10,8}}),
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

      port_a.x=(massFraction/port_a.molarWeight) / AmountOfSolutionInOneKilogram;

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

      Interfaces.SubstanceDefinitionPort port_a
        "The substance with prescribed partial pressure"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      replaceable package substance = Interfaces.BaseSubstance    constrainedby
        Interfaces.BaseSubstance
        "Substance definition: Molar Weight, Enthalpy, Gibbs energy,..."
        annotation (choicesAllMatching = true);

      parameter Boolean usePartialPressureInput = false
        "=true, if fixed partial pressure is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Pressure PartialPressure=0
        "Fixed partial pressure if usePartialPressureInput=false"
        annotation (HideResult=true, Dialog(enable=not usePartialPressureInput));

      parameter Modelica.SIunits.Pressure TotalPressure=101325
        "Total pressure of the whole gaseous solution";

      parameter Modelica.SIunits.Temperature temperature=298.15 "Temperature";
      parameter Modelica.SIunits.MoleFraction moleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential electricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput partialPressure(start=
            PartialPressure, final unit="Pa")=p if usePartialPressureInput
        "Partial pressure of gas = total pressure * gas fraction"
        annotation (HideResult=true,Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.Pressure p "Current partial pressure";

      Modelica.SIunits.AmountOfSubstance amountOfSubstancesInOneLiter =  TotalPressure*0.001/Modelica.Constants.R*temperature
        "Amount of substance in one liter";

       Real activityCoefficient(final unit="1")=
          substance.activityCoefficient(moleFractionBasedIonicStrength,TotalPressure)
        "Activity coefficient of the substance";
    equation
      if not usePartialPressureInput then
        p=PartialPressure;
      end if;

      port_a.x = p / TotalPressure;

      //define the substance at the port
      port_a.activityCoefficient = activityCoefficient;
      port_a.molarWeight = substance.MolarWeight;
      port_a.z = substance.z;

      port_a.molarEnthalpy = substance.molarEnthalpy(electricPotential,TotalPressure);
      port_a.molarEntropy = substance.molarEntropy(port_a.u,temperature,electricPotential,TotalPressure);
      port_a.u0 = substance.u0(temperature,TotalPressure);
      port_a.uPure = substance.uPure(temperature,electricPotential,TotalPressure);
      port_a.molarVolume = substance.molarVolume(TotalPressure,temperature,electricPotential,moleFractionBasedIonicStrength);

      //solution properties at the port
      port_a.temperature = temperature;
      port_a.electricPotential = electricPotential;
      port_a.amountOfSolution = amountOfSubstancesInOneLiter; //particles in one liter of ideal gas

      //electro-chemical potential of the substance in the solution
      port_a.u = substance.u0(temperature,TotalPressure) + Modelica.Constants.R*temperature*log(activityCoefficient* port_a.x)  + substance.z*Modelica.Constants.F*electricPotential;

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

      replaceable package substance = Interfaces.BaseSubstance    constrainedby
        Interfaces.BaseSubstance
        "Substance definition: Molar Weight, Enthalpy, Gibbs energy,..."
        annotation (choicesAllMatching = true);

      parameter Modelica.SIunits.Temperature temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction moleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential electricPotential=0
        "Electric potential";

      Interfaces.SubstanceDefinitionPort port_a
        "constant concentration with any possible flow"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      Real activityCoefficient(final unit="1")=
          substance.activityCoefficient(moleFractionBasedIonicStrength,pressure)
        "Activity coefficient of the substance";

    equation
      port_a.x = 1;

      //define the substance at the port
      port_a.activityCoefficient = activityCoefficient;
      port_a.molarWeight = substance.MolarWeight;
      port_a.z = substance.z;

      port_a.molarEnthalpy = substance.molarEnthalpy(electricPotential,pressure);
      port_a.molarEntropy = substance.molarEntropy(port_a.u,temperature,electricPotential,pressure);
      port_a.u0 = substance.u0(temperature,pressure);
      port_a.uPure = substance.uPure(temperature,electricPotential,pressure);
      port_a.molarVolume = substance.molarVolume(pressure,temperature,electricPotential,moleFractionBasedIonicStrength);

      //solution properties at the port
      port_a.temperature = temperature;
      port_a.electricPotential = electricPotential;
      port_a.amountOfSolution = 1; //particles in one liter of ideal gas

      //electro-chemical potential of the substance in the solution
      port_a.u = port_a.u0 + Modelica.Constants.R*temperature*log(activityCoefficient* port_a.x)  + substance.z*Modelica.Constants.F*electricPotential;

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

      replaceable package substance = Interfaces.BaseSubstance    constrainedby
        Interfaces.BaseSubstance
        "Substance definition: Molar Weight, Enthalpy, Gibbs energy,..."
        annotation (choicesAllMatching = true);

      parameter Modelica.SIunits.Temperature temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction moleFractionBasedIonicStrength=0
        "Ionic strength";
      Modelica.SIunits.ElectricPotential electricPotential "Electric potential";

      Interfaces.SubstanceDefinitionPort port_a
        "constant concentration with any possible flow"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      Real activityCoefficient(final unit="1")=
          substance.activityCoefficient(moleFractionBasedIonicStrength,pressure)
        "Activity coefficient of the substance";

    public
      Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
          Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
              extent={{-110,-10},{-90,10}})));
    equation
      //electric
      pin.v = electricPotential;
      pin.i = -Modelica.Constants.F*port_a.q;

      //pure substance
      port_a.x = 1;

      //define the substance at the port
      port_a.activityCoefficient = activityCoefficient;
      port_a.molarWeight = substance.MolarWeight;
      port_a.z = substance.z;

      port_a.molarEnthalpy = substance.molarEnthalpy(electricPotential,pressure);
      port_a.molarEntropy = substance.molarEntropy(port_a.u,temperature,electricPotential,pressure);
      port_a.u0 = substance.u0(temperature,pressure);
      port_a.uPure = substance.uPure(temperature,electricPotential,pressure);
      port_a.molarVolume = substance.molarVolume(pressure,temperature,electricPotential,moleFractionBasedIonicStrength);

      //solution properties at the port
      port_a.temperature = temperature;
      port_a.electricPotential = electricPotential;
      port_a.amountOfSolution = 1; //particles in one liter of ideal gas

      //electro-chemical potential of the substance in the solution
      port_a.u = port_a.u0 + Modelica.Constants.R*temperature*log(activityCoefficient* port_a.x)  + substance.z*Modelica.Constants.F*electricPotential;

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

       parameter Real Molality(final unit="mol/kg") = 1e-8
        "Fixed molality of the substance if useMolalityInput=false"
        annotation (HideResult=true, Dialog(enable=not useMolalityInput));

       parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionPer1kgOfSolvent = 55.508
        "Amount of all particles in the solution per one kilogram of solvent";

      replaceable package substance = Interfaces.BaseSubstance    constrainedby
        Interfaces.BaseSubstance
        "Substance definition: Molar Weight, Enthalpy, Gibbs energy,..."
        annotation (choicesAllMatching = true);

        parameter Boolean useMolalityInput = false
        "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Temperature temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction moleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential electricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput molalityInput(start=Molality,final unit="mol/kg")=n/KG if
           useMolalityInput
        annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.AmountOfSubstance n "Current amount of the substance";

      Interfaces.SubstanceDefinitionPort port_a
        "constant concentration with any possible flow"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      Real activityCoefficient(final unit="1")=
          substance.activityCoefficient(moleFractionBasedIonicStrength,pressure)
        "Activity coefficient of the substance";
    protected
      constant Modelica.SIunits.Mass KG=1;
    equation
       if not useMolalityInput then
         n=Molality*KG;
       end if;

      port_a.x = n/AmountOfSolutionPer1kgOfSolvent;

      //define the substance at the port
      port_a.activityCoefficient = activityCoefficient;
      port_a.molarWeight = substance.MolarWeight;
      port_a.z = substance.z;

      port_a.molarEnthalpy = substance.molarEnthalpy(electricPotential,pressure);
      port_a.molarEntropy = substance.molarEntropy(port_a.u,temperature,electricPotential,pressure);
      port_a.u0 = substance.u0(temperature,pressure);
      port_a.uPure = substance.uPure(temperature,electricPotential,pressure);
      port_a.molarVolume = substance.molarVolume(pressure,temperature,electricPotential,moleFractionBasedIonicStrength);

      //solution properties at the port
      port_a.temperature = temperature;
      port_a.electricPotential = electricPotential;
      port_a.amountOfSolution = AmountOfSolutionPer1kgOfSolvent; //particles in one liter of ideal gas

      //electro-chemical potential of the substance in the solution
      port_a.u = substance.u0(temperature,pressure) + Modelica.Constants.R*temperature*log(activityCoefficient* port_a.x)  + substance.z*Modelica.Constants.F*electricPotential;

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

       parameter Real MolarConcentration(final unit="mol/m3", displayUnit="mol/l") = 1e-8
        "Fixed molarity of the substance if useMolarityInput=false"
        annotation (HideResult=true, Dialog(enable=not useMolarityInput));

       parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionInOneLiter = 55.508
        "Amount of all particles in the solution one liter of solvent";

       replaceable package substance = Interfaces.BaseSubstance    constrainedby
        Interfaces.BaseSubstance
        "Substance definition: Molar Weight, Enthalpy, Gibbs energy,..."
        annotation (choicesAllMatching = true);

        parameter Boolean useMolarityInput = false
        "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

       parameter Modelica.SIunits.Temperature temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction moleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential electricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput molarConcentrationInput(start=Molality,final unit="mol/m3", displayUnit="mol/l")=n/L if
           useMolarityInput
        annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.AmountOfSubstance n "Current amount of the substance";

      Interfaces.SubstanceDefinitionPort port_a
        "constant concentration with any possible flow"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

       Real activityCoefficient(final unit="1")=
          substance.activityCoefficient(moleFractionBasedIonicStrength)
        "Activity coefficient of the substance";
    protected
      constant Modelica.SIunits.Volume L=0.001;
    equation
       if not useMolarityInput then
         n=MolarConcentration*L;
       end if;

      port_a.x = n/AmountOfSolutionInOneLiter;

       //define the substance at the port
      port_a.activityCoefficient = activityCoefficient;
      port_a.molarWeight = substance.MolarWeight;
      port_a.z = substance.z;

      port_a.molarEnthalpy = substance.molarEnthalpy(electricPotential,pressure);
      port_a.molarEntropy = substance.molarEntropy(port_a.u,temperature,electricPotential,pressure);
      port_a.u0 = substance.u0(temperature,pressure);
      port_a.uPure = substance.uPure(temperature,electricPotential,pressure);
      port_a.molarVolume = substance.molarVolume(pressure,temperature,electricPotential,moleFractionBasedIonicStrength);

      //solution properties at the port
      port_a.temperature = temperature;
      port_a.electricPotential = electricPotential;
      port_a.amountOfSolution = AmountOfSolutionInOneLiter; //particles in one liter of ideal gas

      //electro-chemical potential of the substance in the solution
      port_a.u = substance.u0(temperature,pressure) + Modelica.Constants.R*temperature*log(activityCoefficient* port_a.x)  + substance.z*Modelica.Constants.F*electricPotential;

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

       parameter Modelica.SIunits.MoleFraction MoleFraction = 1e-8
        "Fixed mole fraction of the substance if useMoleFractionInput=false"
        annotation (HideResult=true, Dialog(enable=not useMoleFractionInput));

       parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution = 55.508
        "Amount of all reacting particles in the solution";

      replaceable package substance = Interfaces.BaseSubstance    constrainedby
        Interfaces.BaseSubstance
        "Substance definition: Molar Weight, Enthalpy, Gibbs energy,..."
        annotation (choicesAllMatching = true);

        parameter Boolean useMoleFractionInput = false
        "Is mole fraction of the substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Temperature temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction moleFractionBasedIonicStrength=0
        "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential electricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput moleFractionInput(
        final unit="mol/mol",
        start=MoleFraction)=x if
           useMoleFractionInput annotation (HideResult=true, Placement(transformation(
              extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.MoleFraction x "Current mole fraction of the substance";

      Interfaces.SubstanceDefinitionPort port_a
        "constant concentration with any possible flow"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      Real activityCoefficient(final unit="1")=
          substance.activityCoefficient(moleFractionBasedIonicStrength,pressure)
        "Activity coefficient of the substance";

    equation
       if not useMoleFractionInput then
         x=MoleFraction;
       end if;

      port_a.x = x;

      //define the substance at the port
      port_a.activityCoefficient = activityCoefficient;
      port_a.molarWeight = substance.MolarWeight;
      port_a.z = substance.z;

      port_a.molarEnthalpy = substance.molarEnthalpy(electricPotential,pressure);
      port_a.molarEntropy = substance.molarEntropy(port_a.u,temperature,electricPotential,pressure);
      port_a.u0 = substance.u0(temperature,pressure);
      port_a.uPure = substance.uPure(temperature,electricPotential,pressure);
      port_a.molarVolume = substance.molarVolume(pressure,temperature,electricPotential,moleFractionBasedIonicStrength);

      //solution properties at the port
      port_a.temperature = temperature;
      port_a.electricPotential = electricPotential;
      port_a.amountOfSolution = AmountOfSolution; //particles in one liter of ideal gas

      //electro-chemical potential of the substance in the solution
      port_a.u = substance.u0(temperature,pressure) + Modelica.Constants.R*temperature*log(activityCoefficient* port_a.x)  + substance.z*Modelica.Constants.F*electricPotential;

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

      parameter Modelica.SIunits.VolumeFlowRate Clearance=0
        "Physiological clearance of the substance if useSolutionFlowInput=false"
        annotation (HideResult=true, Dialog(enable=not useSolutionFlowInput));

      parameter Real K(unit="1")=1
        "Coefficient such that Clearance = K*solutionFlow";

      Interfaces.SubstanceUsePort port_a "The substance to be cleared."
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

      Modelica.SIunits.MolarFlowRate molarClearance "Current molar clearance";

    protected
      constant Modelica.SIunits.Volume OneLiter=0.001 "One liter";

    equation
      molarClearance = q*K;

      port_a.q = molarClearance * port_a.x;

      assert(molarClearance>=-Modelica.Constants.eps, "Clearance can not be negative!");

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                            graphics={
            Rectangle(
              extent={{-100,-50},{100,50}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-80,25},{80,0},{-80,-25},{-80,25}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,-100},{150,-60}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-100,-30},{100,-50}},
              lineColor={0,0,0},
              textString="K=%K")}),        Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Clearance;

    model Degradation "Degradation of the substance"

      Interfaces.SubstanceUsePort port_a "Degraded substance"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

      parameter Physiolibrary.Types.Time HalfTime
        "Degradation half time. The time after which will remain half of initial concentration in the defined volume when no other generation, clearence and degradation exist.";

    equation
      port_a.q = (Modelica.Math.log(2)/HalfTime)*port_a.x*port_a.amountOfSolution;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
                            graphics={
            Rectangle(
              extent={{-100,-50},{100,50}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-80,26},{62,0},{-80,-26},{-80,26}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,-100},{150,-60}},
              textString="%name",
              lineColor={0,0,255}),
            Text(
              extent={{-100,-30},{100,-50}},
              lineColor={0,0,0},
              textString="half-time=%HalfTime s"),
            Polygon(
              points={{-68,24},{-68,-24},{-58,-22},{-58,22},{-68,24}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-46,20},{-46,-20},{-36,-18},{-36,18},{-46,20}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-24,16},{-24,-16},{-14,-14},{-14,14},{-24,16}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-2,12},{-2,-12},{8,-10},{8,10},{-2,12}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{20,8},{20,-8},{30,-6},{30,6},{20,8}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{40,4},{40,-4},{50,-2},{50,2},{40,4}},
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
  end Sources;

  package Interfaces
    package BaseSubstance "Base substance"

     constant Modelica.SIunits.MolarMass MolarWeight(displayUnit="kDa")=0
        "Molar weight of the substance in kg/mol or kDa";

     constant Modelica.SIunits.ChargeNumberOfIon z=0
        "Charge number of the substance (e.g. 0..uncharged, -1..electron, +2..Ca^2+)";

     constant Modelica.SIunits.MolarEnergy DfH(displayUnit="kJ/mol")=0
        "Enthalpy of formation of the substance in the selected state";
     constant Modelica.SIunits.MolarEnergy DfG_25degC(displayUnit="kJ/mol")=0
        "Gibbs enerfy of formation at 25°C of the substance in the selected state";

     constant String References[:]={""}
        "References of these thermodynamical values";

     replaceable function activityCoefficient
        "Return activity coefficient of the substance in the solution"
        extends Modelica.Icons.Function;
        input Modelica.SIunits.MoleFraction I
          "Ionic strengh (mole fraction based)";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        output Real gamma "Activity Coefficient";
     algorithm
         gamma := 1;
     end activityCoefficient;

     replaceable function molarEnthalpy "Molar enthalpy of the substance"
        extends Modelica.Icons.Function;
        input Modelica.SIunits.ElectricPotential v
          "Electric potential of the substance";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        output Modelica.SIunits.MolarEnthalpy molarEnthalpy "Molar enthalpy";
     algorithm
         molarEnthalpy := DfH + Modelica.Constants.F*z*v;
     end molarEnthalpy;

     replaceable function molarEntropy "Molar entropy of the substance"
        extends Modelica.Icons.Function;
        input Modelica.SIunits.ChemicalPotential u
          "Electro-chemical potential of the substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        output Modelica.SIunits.MolarEntropy molarEntropy "Molar entropy";
     algorithm
         molarEntropy :=  (u - molarEnthalpy(v,p))/T;
     end molarEntropy;

     replaceable function u0
        "Chemical part of electro-chemical potential of the pure substance"
        extends Modelica.Icons.Function;
        input Modelica.SIunits.Temperature T "Temperature";
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        output Modelica.SIunits.ChemicalPotential u0 "Base chemical potential";
     algorithm
         u0 := DfH - T*((DfH-DfG_25degC)/298.15);
     end u0;

     replaceable function uPure
        "Electro-chemical potential of the pure substance"
        extends Modelica.Icons.Function;
        input Modelica.SIunits.Temperature T "Temperature";
        input Modelica.SIunits.ElectricPotential v=0
          "Electric potential of the substance";
        input Modelica.SIunits.Pressure p=101325 "Pressure";

        output Modelica.SIunits.ChemicalPotential uPure
          "Base electro-chemical potential";
     algorithm
         uPure := u0(T,p) + Modelica.Constants.F*z*v;
     end uPure;

     //not needed for isobaric processes:
     replaceable function molarVolume "Molar volume of the substance"
        extends Modelica.Icons.Function;
        input Modelica.SIunits.Pressure p=101325 "Pressure";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
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
    end BaseSubstance;
    extends Modelica.Icons.InterfacesPackage;

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

      //substance properties (expressed from sunstance definition and current state of the solution)
      output Modelica.SIunits.MolarMass molarWeight(displayUnit="kDa")
        "Molar weight of the substance in kg/mol or kDa";

      output Modelica.SIunits.ChargeNumberOfIon z
        "Charge number of the substance (e.g. 0..uncharged, -1..electron, +2..Ca^2+)";

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
      output Modelica.SIunits.ElectricPotential electricPotential
        "Total electric potential of the solution";

      output Modelica.SIunits.Temperature temperature(displayUnit="degC")
        "Temperature of the solution";

      output Modelica.SIunits.AmountOfSubstance amountOfSolution
        "Total amount of all particles in the solution";

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

      //substance properties (read only; user does not need to calculate or set they again)
      input Modelica.SIunits.MolarMass molarWeight(displayUnit="kDa")
        "Molar weight of the substance in kg/mol or kDa";

      input Modelica.SIunits.ChargeNumberOfIon z
        "Charge number of the substance (e.g. 0..uncharged, -1..electron, +2..Ca^2+)";

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
      input Modelica.SIunits.ElectricPotential electricPotential
        "Total electric potential of the solution";

      input Modelica.SIunits.Temperature temperature(displayUnit="degC")
        "Temperature of the solution";

      input Modelica.SIunits.AmountOfSubstance amountOfSolution
        "Total amount of all particles in the solution";

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

    partial model OnePortParallel
      "Partial molar flow beween two substance definitions"

      SubstanceUsePort port_a
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      SubstanceUsePort port_b
        annotation (Placement(transformation(extent={{90,-10},{110,10}}),
            iconTransformation(extent={{90,-10},{110,10}})));
    equation
      port_a.q + port_b.q = 0;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end OnePortParallel;

    partial model OnePortSerial
      "Partial transfer of substance from substance definition component to another transfer component (such as MolarFlowSensor)"

      SubstanceUsePort port_a
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      SubstanceDefinitionPort port_b
        annotation (Placement(transformation(extent={{90,-10},{110,10}}),
            iconTransformation(extent={{90,-10},{110,10}})));
    equation
      port_a.q + port_b.q = 0;

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
      port_b.electricPotential = port_a.electricPotential;
      port_b.amountOfSolution = port_a.amountOfSolution;

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
</html>"));
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
      solution.n = amountOfSolution;

      //heat
      der(solution.S) = solution.dS;
      der(freeEnthalpy) = solution.dH + heatFromEnvironment;
      der(solution.G) = solution.dG;
      solution.T = (freeEnthalpy-solution.G)/solution.S;

      //electric
      der(charge) = solution.i;

      //ionic strength (mole fraction based)
      der(solution.I*amountOfSolution) = solution.dI;

      //isobaric
      der(volume) = solution.dV;

                                                                                                          annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false, initialScale=1, extent={{-100,-100},{100,100}}),
            graphics={Text(
              extent={{-88,-92},{92,-100}},
              lineColor={0,0,255},
              textString="%name",
              horizontalAlignment=TextAlignment.Left),
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
