within Chemical;
package Examples "Examples that demonstrate usage of chemical library"
extends Modelica.Icons.ExamplesPackage;

  model SimpleReaction
    "The simple chemical reaction A<->B with equilibrium B/A = 2"
     extends Modelica.Icons.Example;

    constant Real K = 2 "Dissociation constant of the reaction";

    constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
    constant Real R = Modelica.Constants.R "Gas constant";

    Chemical.Components.Solution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    Chemical.Components.Substance A(
      substanceData(MolarWeight=1),
      use_mass_start=false,
      amountOfSubstance_start=0.9)
      annotation (Placement(transformation(extent={{-52,-8},{-32,12}})));

    Chemical.Components.Reaction reaction(nS=1, nP=1)
      annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
    Chemical.Components.Substance B(
      substanceData(DfG=-R*T_25degC*log(K), MolarWeight=1),
      use_mass_start=false,
      amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{62,-8},{42,12}})));

  equation
    connect(A.solution, solution.solution) annotation (Line(
        points={{-48,-8},{-48,-92},{60,-92},{60,-98}},
        color={127,127,0}));
    connect(B.solution, solution.solution) annotation (Line(points={{58,-8},{
          58,-92},{60,-92},{60,-98}},  color={127,127,0}));
    connect(A.port_a, reaction.substrates[1]) annotation (Line(points={{-32,2},
            {-22,2},{-22,2},{-10,2}}, color={158,66,200}));
    connect(reaction.products[1], B.port_a)
      annotation (Line(points={{10,2},{42,2}}, color={158,66,200}));
    annotation (Documentation(revisions="<html>
<p><i>2015-2018</i></p>
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

    constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
    constant Real R = Modelica.Constants.R "Gas constant";

    Chemical.Components.Solution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    Chemical.Components.Substance A(use_mass_start=false, amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
    Chemical.Components.Reaction reaction(nS=2, nP=1)
      annotation (Placement(transformation(extent={{4,-8},{24,12}})));
    Chemical.Components.Substance B(use_mass_start=false, amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
    Chemical.Components.Substance C(
      substanceData(DfG=-R*T_25degC*log(Kx)),
      use_mass_start=false,
      amountOfSubstance_start=0.1)
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
        points={{-14,-14},{-10,-14},{-10,4},{4,4}},
        color={158,66,200},
        thickness=1));
    connect(A.port_a, reaction.substrates[2]) annotation (Line(
        points={{-14,12},{-10,12},{-10,0},{4,0}},
        color={158,66,200},
        thickness=1));
    annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
      experiment(StopTime=0.0001, __Dymola_Algorithm="Dassl"));
  end SimpleReaction2;

  model HeatingOfWater "Heating of 1 kg water"
    extends Modelica.Icons.Example;

    Chemical.Components.Solution solution(
      useMechanicPorts=true,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
      annotation (Placement(transformation(extent={{-86,-72},{-66,-52}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
    Components.Substance liquidWater(mass_start=1, substanceData=
          Chemical.Substances.Water_liquid())
      annotation (Placement(transformation(extent={{22,-28},{42,-8}})));
    inner Modelica.Fluid.System system(T_ambient=298.15)
      annotation (Placement(transformation(extent={{60,50},{80,70}})));
  equation
    connect(fixedHeatFlow.port, solution.heatPort) annotation (Line(
        points={{-66,-62},{-60,-62},{-60,-102}},
        color={191,0,0}));
  connect(fixed1.flange, solution.bottom) annotation (Line(
      points={{-18,-84},{0,-84},{0,-102}},
      color={0,127,0}));
    connect(solution.solution, liquidWater.solution) annotation (Line(points={{
            60,-98},{26,-98},{26,-28}}, color={127,127,0}));
    annotation (experiment(StopTime=1),
    Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Heating of solution by one degree, using standard HeatPort from Modelica Standard Library.</p>
<p>Observe Solution.T (or H2O.Solution.T) for temperature change.</p>
</html>"));
  end HeatingOfWater;

  model HeatingOfAlcohol "Heating of 50% ethanol"
    extends Modelica.Icons.Example;

    Chemical.Components.Solution solution(
      useMechanicPorts=true,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
      annotation (Placement(transformation(extent={{-86,-76},{-66,-56}})));
    Chemical.Components.Substance Ethanol(
      redeclare package stateOfMatter =
          Chemical.Interfaces.Incompressible,
      substanceData=Chemical.Substances.Ethanol_liquid(),
      mass_start=(55.508/2)*0.04607)
      annotation (Placement(transformation(extent={{18,-8},{38,12}})));

    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
    Components.Substance liquidWater(mass_start=1/2, substanceData=
          Chemical.Substances.Water_liquid())
      annotation (Placement(transformation(extent={{-50,-8},{-30,12}})));
    inner Modelica.Fluid.System system(T_ambient=298.15)
      annotation (Placement(transformation(extent={{62,42},{82,62}})));
  equation
    connect(fixedHeatFlow.port, solution.heatPort) annotation (Line(
        points={{-66,-66},{-60,-66},{-60,-102}},
        color={191,0,0}));
  connect(solution.solution, Ethanol.solution) annotation (Line(
      points={{60,-98},{60,-34},{22,-34},{22,-8}},
      color={127,127,0}));
  connect(solution.bottom, fixed1.flange) annotation (Line(
      points={{0,-102},{0,-84},{-18,-84}},
      color={0,127,0}));
    connect(solution.solution, liquidWater.solution) annotation (Line(points={{
            60,-98},{60,-34},{-46,-34},{-46,-8}}, color={127,127,0}));
    annotation (experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Heating of solution of water and ethanol, using standard HeatPort from Modelica Standard Library.</p>
<p>Observe Solution.T (or H2O.Solution.T or Ethanol.Solution.T) for temperature change. Note, that we can heat all substances in solution at once and the results would differ from HeatingOfWater.</p>
</html>"));
  end HeatingOfAlcohol;

  model ExothermicReaction
    "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"

     extends Modelica.Icons.Example;

    parameter Modelica.Units.SI.MolarEnergy ReactionEnthalpy=-55000;

    Chemical.Components.Solution thermal_isolated_solution(
      useMechanicPorts=true,
      ConstantTemperature=false)
      annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));
    Chemical.Components.Substance A(use_mass_start=false, amountOfSubstance_start=0.9)
      annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
    Chemical.Components.Reaction reaction(nS=1, nP=1)
      annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
    Chemical.Components.Substance B(
      substanceData(DfH=ReactionEnthalpy),
      use_mass_start=false,
      amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{40,-60},{20,-40}})));

    Chemical.Components.Solution solution_at_constant_temperature(
      useMechanicPorts=true,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{-100,0},{98,94}})));
    Chemical.Components.Substance A1(use_mass_start=false, amountOfSubstance_start=0.9)
      annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
    Chemical.Components.Reaction reaction1(nS=1, nP=1)
      annotation (Placement(transformation(extent={{-8,40},{12,60}})));
    Chemical.Components.Substance B1(
      substanceData(DfH=ReactionEnthalpy),
      use_mass_start=false,
      amountOfSubstance_start=0.1)
      annotation (Placement(transformation(extent={{40,40},{20,60}})));

    //  Modelica.SIunits.HeatFlowRate q
    //    "Heat flow to environment to reach constant temperature";
    Modelica.Units.SI.Temperature t
      "Temperature if the solution is ideally thermal isolated from environment";
    Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
        mass_start=1)
      annotation (Placement(transformation(extent={{20,4},{40,24}})));
    Components.Substance H2O1(substanceData=Chemical.Substances.Water_liquid(),
        mass_start=1)
      annotation (Placement(transformation(extent={{20,-94},{40,-74}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-28,4},{-8,24}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed2
      annotation (Placement(transformation(extent={{-26,-96},{-6,-76}})));
    inner Modelica.Fluid.System system(T_ambient=298.15)
      annotation (Placement(transformation(extent={{56,64},{76,84}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=
          298.15)
      annotation (Placement(transformation(extent={{-88,26},{-68,46}})));
  equation
    //  q = fixedTemperature.port.Q_flow;
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
    connect(solution_at_constant_temperature.heatPort, fixedTemperature.port)
      annotation (Line(points={{-60.4,-0.94},{-60.4,36},{-68,36}}, color={191,0,
            0}));
    annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Demonstration of exotermic reaction with perfect cooling (i.e. connected fixed temperature to the HeatPort) and thermally insulated (HetPort unconnected). See solution_(...).T</p>
</html>"),
      experiment(StopTime=0.001),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}})));
  end ExothermicReaction;

  model HydrogenCombustion "Hydrogen combustion in piston"
    extends Modelica.Icons.Example;

    parameter Modelica.Units.SI.Volume V=0.001 "Initial volume";
   // parameter Modelica.SIunits.Pressure p=100000 "Initial pressure";
    parameter Modelica.Units.SI.Temperature T=298.15
      "Initial temperature";

    parameter Modelica.Units.SI.Area A=0.01 "Cross area of cylinder";

    //p*V=n*R*T
    // parameter Modelica.SIunits.AmountOfSubstance n=p*V/(Modelica.Constants.R*T)
    //   "Initial amount of substances in sulution";
    Chemical.Components.Solution idealGas(
      SurfaceArea=A,
      useMechanicPorts=true,
      useThermalPort=true,
      redeclare package stateOfMatter = Interfaces.IdealGas)
      annotation (Placement(transformation(extent={{-108,-50},{-8,50}})));
                     // AmbientPressure=p)
    //  volume_start=V,
    Chemical.Components.Substance H2_gas(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Hydrogen_gas(),
      use_mass_start=false,
      amountOfSubstance_start=0.026)
      annotation (Placement(transformation(extent={{-98,-26},{-78,-6}})));
    Chemical.Components.Substance O2_gas(
      substanceData=Chemical.Substances.Oxygen_gas(),
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      use_mass_start=false,
      amountOfSubstance_start=0.013)
      annotation (Placement(transformation(extent={{-98,10},{-78,30}})));
    Chemical.Components.Substance H2O_gas(
      substanceData=Chemical.Substances.Water_gas(),
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      use_mass_start=false,
      amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-14,-8},{-34,12}})));
    Chemical.Components.Reaction reaction(
      nS=2,
      s={2,1},
      p={2},
      kE=4e-05,
      nP=1) annotation (Placement(transformation(extent={{-68,-8},{-48,12}})));
    Modelica.Mechanics.Translational.Components.Spring spring(c=1e6) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-58,64})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=2)
      annotation (Placement(transformation(extent={{-98,-80},{-78,-60}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature coolerTemperature(T=298.15)
      annotation (Placement(transformation(extent={{-18,-80},{-38,-60}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-58,78})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-68,-66},{-48,-46}})));
    Components.Solution          idealGas1(
      SurfaceArea=A,
      useMechanicPorts=true,
      useThermalPort=true,
      redeclare package stateOfMatter = Interfaces.IdealGas)
      annotation (Placement(transformation(extent={{18,-52},{118,48}})));
    Components.Substance H2_gas1(
      redeclare package stateOfMatter = Interfaces.IdealGasMSL,
      substanceData(data=Modelica.Media.IdealGases.Common.SingleGasesData.H2),
      use_mass_start=false,
      amountOfSubstance_start=0.026)
      annotation (Placement(transformation(extent={{28,-28},{48,-8}})));
    Components.Substance O2_gas1(
      substanceData(data=Modelica.Media.IdealGases.Common.SingleGasesData.O2),
      redeclare package stateOfMatter = Interfaces.IdealGasMSL,
      use_mass_start=false,
      amountOfSubstance_start=0.013)
      annotation (Placement(transformation(extent={{28,8},{48,28}})));
    Components.Substance H2O_gas1(
      substanceData(data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O),
      redeclare package stateOfMatter = Interfaces.IdealGasMSL,
      use_mass_start=false,
      amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{112,-10},{92,10}})));

    Components.Reaction          reaction1(
      nS=2,
      s={2,1},
      p={2},
      kE=4e-05,
      nP=1) annotation (Placement(transformation(extent={{58,-10},{78,10}})));
    Modelica.Mechanics.Translational.Components.Spring spring1(c=1e6)
                                                                     annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={68,62})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(G=2)
      annotation (Placement(transformation(extent={{28,-82},{48,-62}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature coolerTemperature1(T=298.15)
      annotation (Placement(transformation(extent={{108,-82},{88,-62}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed2
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={68,76})));
    Modelica.Mechanics.Translational.Components.Fixed fixed3
      annotation (Placement(transformation(extent={{58,-68},{78,-48}})));
  equation
  connect(H2_gas.port_a, reaction.substrates[1]) annotation (Line(
      points={{-78,-16},{-74,-16},{-74,4},{-68,4}},
      color={158,66,200},
      thickness=1));
  connect(O2_gas.port_a, reaction.substrates[2]) annotation (Line(
      points={{-78,20},{-74,20},{-74,0},{-68,0}},
      color={158,66,200},
      thickness=1));
  connect(H2_gas.solution, idealGas.solution) annotation (Line(
      points={{-94,-26},{-28,-26},{-28,-49}},
      color={127,127,0}));
  connect(O2_gas.solution, idealGas.solution) annotation (Line(
      points={{-94,10},{-102,10},{-102,-26},{-28,-26},{-28,-49}},
      color={127,127,0}));
  connect(H2O_gas.solution, idealGas.solution) annotation (Line(
      points={{-18,-8},{-18,-26},{-28,-26},{-28,-49}},
      color={127,127,0}));
    connect(idealGas.surfaceFlange, spring.flange_a) annotation (Line(
        points={{-58,50},{-58,54}},
        color={0,127,0}));
    connect(idealGas.heatPort, thermalConductor.port_a) annotation (Line(
        points={{-88,-51},{-88,-56},{-106,-56},{-106,-70},{-98,-70}},
        color={191,0,0}));
    connect(thermalConductor.port_b, coolerTemperature.port) annotation (Line(
        points={{-78,-70},{-38,-70}},
        color={191,0,0}));
    connect(fixed.flange, spring.flange_b) annotation (Line(
        points={{-58,78},{-58,74}},
        color={0,127,0}));
  connect(idealGas.bottom, fixed1.flange) annotation (Line(
      points={{-58,-51},{-58,-56}},
      color={0,127,0}));
    connect(reaction.products[1], H2O_gas.port_a) annotation (Line(points={{-48,2},
            {-34,2}},                     color={158,66,200}));
    connect(H2_gas1.port_a, reaction1.substrates[1]) annotation (Line(
        points={{48,-18},{52,-18},{52,2},{58,2}},
        color={158,66,200},
        thickness=1));
    connect(O2_gas1.port_a, reaction1.substrates[2]) annotation (Line(
        points={{48,18},{52,18},{52,-2},{58,-2}},
        color={158,66,200},
        thickness=1));
    connect(H2_gas1.solution, idealGas1.solution) annotation (Line(points={{32,-28},
            {98,-28},{98,-51}},   color={127,127,0}));
    connect(O2_gas1.solution, idealGas1.solution) annotation (Line(points={{32,8},{
            24,8},{24,-28},{98,-28},{98,-51}},    color={127,127,0}));
    connect(H2O_gas1.solution, idealGas1.solution) annotation (Line(points={{108,-10},
            {108,-28},{98,-28},{98,-51}},   color={127,127,0}));
    connect(idealGas1.surfaceFlange, spring1.flange_a)
      annotation (Line(points={{68,48},{68,52}},   color={0,127,0}));
    connect(idealGas1.heatPort, thermalConductor1.port_a) annotation (Line(points={{38,-53},
            {38,-58},{20,-58},{20,-72},{28,-72}},          color={191,0,0}));
    connect(thermalConductor1.port_b, coolerTemperature1.port)
      annotation (Line(points={{48,-72},{88,-72}},   color={191,0,0}));
    connect(fixed2.flange, spring1.flange_b)
      annotation (Line(points={{68,76},{68,72}},   color={0,127,0}));
    connect(idealGas1.bottom, fixed3.flange)
      annotation (Line(points={{68,-53},{68,-58}},   color={0,127,0}));
    connect(reaction1.products[1], H2O_gas1.port_a)
      annotation (Line(points={{78,0},{92,0}},     color={158,66,200}));
    annotation ( experiment(StopTime=0.39, __Dymola_Algorithm="Dassl"),
                                         Documentation(info="<html>
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
</html>"),
      Diagram(coordinateSystem(extent={{-120,-100},{120,100}})));
  end HydrogenCombustion;

  model WaterVaporization "Evaporation of water"
     extends Modelica.Icons.Example;

    parameter Modelica.Units.SI.Temperature T_start=273.15
      "Initial temperature";

    Chemical.Components.Solution liquid(
      temperature_start=T_start,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{-98,-98},{-6,-8}})));

    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
    Chemical.Components.Solution gas(
      temperature_start=T_start,
      useThermalPort=true,
      redeclare package stateOfMatter = Interfaces.IdealGas)
      annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                  /*volume_start(
        displayUnit="l") = 0.001, */
    Components.Substance H2O_gaseuous(
      redeclare package stateOfMatter = Interfaces.IdealGas,
      substanceData=Chemical.Substances.Water_gas(),
      mass_start=0.000106537)
      annotation (Placement(transformation(extent={{28,50},{8,70}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                         fixedTemperature
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        origin={84,8})));
    Modelica.Blocks.Sources.ContinuousClock clock(offset=1*T_start)
      annotation (Placement(transformation(extent={{62,36},{82,56}})));
    Chemical.Components.Substance otherSubstances(
      redeclare package stateOfMatter = Interfaces.IdealGas,
      substanceData=Chemical.Substances.Oxygen_gas(),
      use_mass_start=false,
      amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{2,28},{22,48}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(
      G=1e6) annotation (Placement(transformation(extent={{48,-8},{68,12}})));
    Components.Substance liquidWater(substanceData=
          Chemical.Substances.Water_liquid(), mass_start=1)
      annotation (Placement(transformation(extent={{-28,-62},{-48,-42}})));
    inner Modelica.Fluid.System system(p_ambient=100000, T_ambient=298.15)
      annotation (Placement(transformation(extent={{54,-48},{74,-28}})));
    Components.GasSolubility          gasSolubility
      annotation (Placement(transformation(extent={{-92,16},{-72,36}})));
    Sensors.PartialPressureSensor pH2O(redeclare package stateOfMatter =
          Interfaces.IdealGas, substanceData=
          Chemical.Substances.Water_gas())
      annotation (Placement(transformation(extent={{-26,76},{-6,96}})));
  equation

    connect(gas.solution, H2O_gaseuous.solution) annotation (Line(
        points={{27.6,6.9},{24,6.9},{24,50}},
        color={127,127,0}));
  connect(fixedTemperature.T, clock.y) annotation (Line(
      points={{96,8},{98,8},{98,46},{83,46}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(gas.solution, otherSubstances.solution) annotation (Line(
      points={{27.6,6.9},{6,6.9},{6,28}},
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
    connect(liquid.solution, liquidWater.solution) annotation (Line(points={{-24.4,
            -97.1},{-24.4,-96.55},{-32,-96.55},{-32,-62}}, color={127,127,0}));
    connect(H2O_gaseuous.port_a,gasSolubility. gas_port) annotation (Line(
        points={{8,60},{-82,60},{-82,36}},
        color={158,66,200},
        thickness=1));
    connect(gasSolubility.liquid_port, liquidWater.port_a) annotation (Line(
          points={{-82,16},{-82,-52},{-48,-52}}, color={158,66,200}));
    connect(pH2O.port_a, H2O_gaseuous.port_a) annotation (Line(points={{
            -6,86},{0,86},{0,60},{8,60}}, color={158,66,200}));
    connect(pH2O.solution, gas.solution) annotation (Line(points={{-22,76},
            {-22,6.9},{27.6,6.9}}, color={127,127,0}));
    annotation (
      experiment(
        StopTime=100,
        Tolerance=1e-08,
        __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>Demonstraiton of water vaporization between two solutions - liquid and gaseous. The temperature is increased in time to illustrate, how the vaporization rate rises in higher temperatures. See liquid.T and liquid.Volume, compared to gas.T and gas.volume.</p>
</html>",
        revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end WaterVaporization;

  model WaterSublimation "Sublimation of water"
     extends Modelica.Icons.Example;

    parameter Modelica.Units.SI.Temperature T_start=273.15 - 50
      "Initial temperature";

    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
    Chemical.Components.Solution gas(
      temperature_start=T_start,
      useThermalPort=true,
      redeclare package stateOfMatter = Interfaces.IdealGas,
      BasePressure=600)
      annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                  /*volume_start(
        displayUnit="l") = 0.001, */
    Chemical.Components.Substance H2O_gaseuous(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Water_gas(),
      use_mass_start=false,
      amountOfSubstance_start=0.001)
      annotation (Placement(transformation(extent={{24,56},{4,76}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                         fixedTemperature
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        origin={84,8})));
    Modelica.Blocks.Sources.ContinuousClock clock(offset=1*T_start)
      annotation (Placement(transformation(extent={{62,36},{82,56}})));
    Chemical.Components.Substance otherSubstances(
      substanceData=Chemical.Substances.Oxygen_gas(),
      redeclare package stateOfMatter = Interfaces.IdealGas,
      use_mass_start=false,
      amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-4,36},{16,56}})));
    Chemical.Components.Solution solid(
      temperature_start=T_start,
      BasePressure=600,
      useThermalPort=true)
      annotation (Placement(transformation(extent={{8,-98},{100,-8}})));
    Chemical.Components.Substance H2O_solid(
      substanceData=Chemical.Substances.Water_IceIh(),
      use_mass_start=false,
      amountOfSubstance_start=55.508) "Solid water"
      annotation (Placement(transformation(extent={{70,-62},{50,-42}})));
    Chemical.Components.GasSolubility gasSolubility1(KC=10)
      annotation (Placement(transformation(extent={{-76,18},{-56,38}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=1e6)
      annotation (Placement(transformation(extent={{48,-8},{68,12}})));
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
    connect(fixedTemperature.port, thermalConductor.port_b) annotation (Line(
        points={{74,8},{72,8},{72,2},{68,2}},
        color={191,0,0}));
    connect(gas.heatPort, thermalConductor.port_a) annotation (Line(
        points={{-27.6,5.1},{-28,5.1},{-28,2},{48,2}},
        color={191,0,0}));
    connect(solid.heatPort, thermalConductor.port_a) annotation (Line(
        points={{26.4,-98.9},{-2,-98.9},{-2,2},{48,2}},
        color={191,0,0}));
    annotation (
      experiment(StopTime=49.7, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>Demonstraiton of water sublimation between two solutions - solid and gaseous. The temperature is increased in time to illustrate, how the sublimation rate rises in higher temperatures. See solid.T and solid.Volume, compared to gas.T and gas.volume. Note, that the liquid phase is omitted here.</p>
</html>",
        revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end WaterSublimation;

  model GasSolubility_NIST "Dissolution of gases in liquids"
     extends Modelica.Icons.Example;

    Chemical.Components.Solution water_solution_25degC(temperature_start=298.15)
      annotation (Placement(transformation(extent={{-160,-78},{-68,12}})));
                                        //(amountOfSolution_start=52.3)
    Chemical.Components.Solution water_solution_37degC(temperature_start=310.15)
      annotation (Placement(transformation(extent={{-52,-80},{42,12}})));
                                     //(amountOfSolution_start=39.7)
    Chemical.Components.GasSolubility CO2_dissolutionP
      annotation (Placement(transformation(extent={{-138,42},{-118,62}})));
    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
    Chemical.Components.Substance CO2_25(
      amountOfSubstance_start(displayUnit="mmol") = 0.001,
      substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
      use_mass_start=false) "Free dissolved CO2 in water at 25 degC"
      annotation (Placement(transformation(extent={{-150,-26},{-130,-6}})));
    Chemical.Components.GasSolubility O2_dissolutionP
      annotation (Placement(transformation(extent={{-100,42},{-80,62}})));

    Chemical.Sources.ExternalIdealGasSubstance O2_g_25(
      substanceData=Chemical.Substances.Oxygen_gas(),
      TotalPressure=system.p_ambient,
      PartialPressure(displayUnit="mmHg") = 12665.626804425,
      Temperature=298.15)
      annotation (Placement(transformation(extent={{-114,74},{-94,94}})));
    Chemical.Components.Substance O2_25(
      amountOfSubstance_start(displayUnit="mmol") = 0.0001,
      substanceData=Chemical.Substances.Oxygen_aqueous(),
      use_mass_start=false) "Free dissolved O2 in water at 25 degC"
      annotation (Placement(transformation(extent={{-114,-26},{-94,-6}})));
    Chemical.Components.GasSolubility CO2_dissolutionE
      annotation (Placement(transformation(extent={{-26,42},{-6,62}})));

    Chemical.Sources.ExternalIdealGasSubstance CO2_g_25(
      substanceData=Chemical.Substances.CarbonDioxide_gas(),
      TotalPressure=system.p_ambient,
      PartialPressure(displayUnit="mmHg") = 5332.8954966,
      Temperature=298.15)
      annotation (Placement(transformation(extent={{-154,74},{-134,94}})));

    Chemical.Components.Substance CO2_37(
      amountOfSubstance_start(displayUnit="mmol") = 0.001,
      substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
      use_mass_start=false) "Free dissolved CO2 in water at 37degC"
      annotation (Placement(transformation(extent={{-42,-34},{-22,-14}})));

    Chemical.Components.GasSolubility O2_dissolutionE_NIST
      annotation (Placement(transformation(extent={{18,42},{38,62}})));
    Chemical.Components.Substance O2_37(
      amountOfSubstance_start(displayUnit="mmol") = 0.0001,
      substanceData=Chemical.Substances.Oxygen_aqueous(),
      use_mass_start=false) "Free dissolved O2 in water at 37degC"
      annotation (Placement(transformation(extent={{-2,-34},{18,-14}})));
    Components.Substance water_25(substanceData=
          Chemical.Substances.Water_liquid(), mass_start=1)
      annotation (Placement(transformation(extent={{-100,-68},{-80,-48}})));
    Components.Substance water_37(substanceData=
          Chemical.Substances.Water_liquid(), mass_start=1)
      annotation (Placement(transformation(extent={{8,-70},{28,-50}})));
    Sources.ExternalIdealGasSubstance CO2_g_37(
      substanceData=Chemical.Substances.CarbonDioxide_gas(),
      TotalPressure=system.p_ambient,
      PartialPressure(displayUnit="mmHg") = 5332.8954966,
      Temperature=310.15)
      annotation (Placement(transformation(extent={{-44,68},{-24,88}})));
    Sources.ExternalIdealGasSubstance O2_g_37(
      substanceData=Chemical.Substances.Oxygen_gas(),
      TotalPressure=system.p_ambient,
      PartialPressure(displayUnit="mmHg") = 12665.626804425,
      Temperature=310.15)
      annotation (Placement(transformation(extent={{-6,68},{14,88}})));
    Components.Solution          water_solution_37degC1(temperature_start=
          273.15)
      annotation (Placement(transformation(extent={{66,-80},{160,12}})));
    Components.GasSolubility          CO2_dissolutionE1
      annotation (Placement(transformation(extent={{92,42},{112,62}})));
    Components.Substance CO2_0(
      amountOfSubstance_start(displayUnit="mmol") = 0.001,
      substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
      use_mass_start=false) "Free dissolved CO2 in water at 0degC"
      annotation (Placement(transformation(extent={{76,-34},{96,-14}})));
    Components.GasSolubility          O2_dissolutionE_NIST1
      annotation (Placement(transformation(extent={{136,42},{156,62}})));
    Components.Substance O2_0(
      amountOfSubstance_start(displayUnit="mmol") = 0.0001,
      substanceData=Chemical.Substances.Oxygen_aqueous(),
      use_mass_start=false) "Free dissolved O2 in water at 0degC"
      annotation (Placement(transformation(extent={{116,-34},{136,-14}})));
    Components.Substance water_0(substanceData=Chemical.Substances.Water_liquid(),
        mass_start=1)
      annotation (Placement(transformation(extent={{126,-70},{146,-50}})));
    Sources.ExternalIdealGasSubstance CO2_g_0(
      substanceData=Chemical.Substances.CarbonDioxide_gas(),
      TotalPressure=system.p_ambient,
      PartialPressure(displayUnit="mmHg") = 5332.8954966,
      Temperature=273.15)
      annotation (Placement(transformation(extent={{74,68},{94,88}})));
    Sources.ExternalIdealGasSubstance O2_g_0(
      substanceData=Chemical.Substances.Oxygen_gas(),
      TotalPressure=system.p_ambient,
      PartialPressure(displayUnit="mmHg") = 12665.626804425,
      Temperature=273.15)
      annotation (Placement(transformation(extent={{112,68},{132,88}})));
    inner Modelica.Fluid.System system(p_ambient=100000)
      annotation (Placement(transformation(extent={{-70,-98},{-50,-78}})));
    Real kH_CO2_25, kH_O2_25;
  equation

    kH_CO2_25 = CO2_25.c / CO2_g_25.x;
    kH_O2_25 = O2_25.c / O2_g_25.x;
  //  kH_CO2_25 = CO2_25.x / CO2_g_25.x;
  //  kH_CO2_25 = CO2_25.x / CO2_g_25.x;


  connect(CO2_g_25.port_a, CO2_dissolutionP.gas_port) annotation (Line(
      points={{-134,84},{-128,84},{-128,62}},
      color={158,66,200},
      thickness=1));
    connect(CO2_dissolutionP.liquid_port, CO2_25.port_a) annotation (Line(
        points={{-128,42},{-128,-16},{-130,-16}},
        color={158,66,200},
        thickness=1));
    connect(CO2_dissolutionE.liquid_port, CO2_37.port_a) annotation (Line(
        points={{-16,42},{-16,-24},{-22,-24}},
        color={158,66,200},
        thickness=1));
  connect(O2_g_25.port_a, O2_dissolutionP.gas_port) annotation (Line(
      points={{-94,84},{-90,84},{-90,62}},
      color={158,66,200},
      thickness=1));
    connect(O2_dissolutionP.liquid_port, O2_25.port_a) annotation (Line(
        points={{-90,42},{-90,-14},{-94,-14},{-94,-16}},
        color={158,66,200},
        thickness=1));
    connect(CO2_25.solution, water_solution_25degC.solution) annotation (Line(
          points={{-146,-26},{-146,-77.1},{-86.4,-77.1}},
                                                        color={127,127,0}));
    connect(O2_25.solution, water_solution_25degC.solution) annotation (Line(
          points={{-110,-26},{-110,-77.1},{-86.4,-77.1}},
                                                        color={127,127,0}));
    connect(CO2_37.solution, water_solution_37degC.solution) annotation (Line(
          points={{-38,-34},{-38,-79.08},{23.2,-79.08}},
                                                       color={127,127,0}));
    connect(O2_dissolutionE_NIST.liquid_port, O2_37.port_a) annotation (Line(
        points={{28,42},{28,-24},{18,-24}},
        color={158,66,200},
        thickness=1));
    connect(O2_37.solution, water_solution_37degC.solution) annotation (Line(
          points={{2,-34},{2,-79.08},{23.2,-79.08}},   color={127,127,0}));
    connect(water_25.solution, water_solution_25degC.solution) annotation (Line(
          points={{-96,-68},{-96,-77.1},{-86.4,-77.1}}, color={127,127,0}));
    connect(water_37.solution, water_solution_37degC.solution) annotation (Line(
          points={{12,-70},{12,-79.08},{23.2,-79.08}}, color={127,127,0}));
    connect(CO2_g_37.port_a, CO2_dissolutionE.gas_port)
      annotation (Line(points={{-24,78},{-16,78},{-16,62}},
                                                         color={158,66,200}));
    connect(O2_g_37.port_a, O2_dissolutionE_NIST.gas_port)
      annotation (Line(points={{14,78},{28,78},{28,62}}, color={158,66,200}));
    connect(CO2_dissolutionE1.liquid_port, CO2_0.port_a) annotation (Line(
        points={{102,42},{102,-24},{96,-24}},
        color={158,66,200},
        thickness=1));
    connect(CO2_0.solution, water_solution_37degC1.solution) annotation (Line(
          points={{80,-34},{80,-79.08},{141.2,-79.08}},   color={127,127,0}));
    connect(O2_dissolutionE_NIST1.liquid_port, O2_0.port_a) annotation (Line(
        points={{146,42},{146,-24},{136,-24}},
        color={158,66,200},
        thickness=1));
    connect(O2_0.solution, water_solution_37degC1.solution) annotation (Line(
          points={{120,-34},{120,-79.08},{141.2,-79.08}}, color={127,127,0}));
    connect(water_0.solution, water_solution_37degC1.solution) annotation (Line(
          points={{130,-70},{130,-79.08},{141.2,-79.08}}, color={127,127,0}));
    connect(CO2_g_0.port_a, CO2_dissolutionE1.gas_port) annotation (Line(points={{94,78},
            {102,78},{102,62}},          color={158,66,200}));
    connect(O2_g_0.port_a, O2_dissolutionE_NIST1.gas_port) annotation (Line(
          points={{132,78},{146,78},{146,62}}, color={158,66,200}));
    annotation (
      experiment(StopTime=1e-005),
      Documentation(info="<html>
<p>Demonstration of CO2 and O2 dissolution in pure water as described by NIST Henry&apos;s law data at 25degC.</p>
<p>Recalculation from Henry&apos;s law constants ( https://webbook.nist.gov/cgi/inchi?ID=C124389&amp;Mask=10#Solubility ):</p>
<p>CO2:</p>
<p>kH = 0.035 mol/(kg.bar) at 25degC</p>
<ul>
<li>pCO2 = 40 mmHg =&gt; dissolved CO2 .. 1.87 mmol/kg at 25degC</li>
</ul>
<p><br>Other temperatures:</p>
<p>NIST constant ... 2400 K ( https://webbook.nist.gov/cgi/inchi?ID=C124389&amp;Mask=10#Solubility )</p>
<p>0degC ... 0.035*exp(2400*(1/273.15 - 1/298.15) = 0.073 mol/(kg.bar)</p>
<ul>
<li>pCO2 = 40 mmHg =&gt; dissolved CO2 .. 3.9 mmol/kg at 0degC</li>
</ul>
<p>37degC ... 0.035*exp(2400*(1/310.15 - 1/298.15) = 0.026 mol/(kg.bar)</p>
<ul>
<li>pCO2 = 40 mmHg =&gt; dissolved CO2 .. 1.4 mmol/kg at 37degC</li>
</ul>
<p><br>O2:</p>
<p>kH = 0.0013 mol/(kg.bar) at 25degC</p>
<ul>
<li>pO2 = 95 mmHg (0.126656 bar) =&gt; dissolved O2 .. 0.165 mmol/kg at 25degC</li>
</ul>
<p><br>Other temperatures:</p>
<p>NIST constant ... 1500 K ( https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&amp;Mask=10 )</p>
<p>0degC ... 0.0013*exp(1500*(1/273.15 - 1/298.15) = 0.0021 mol/(kg.bar)</p>
<ul>
<li>pO2 = 95 mmHg =&gt; dissolved O2 .. 0.26 mmol/kg at 0degC</li>
</ul>
<p>37degC ... 0.0013*exp(1500*(1/273.15 - 1/298.15) = 0.0011 mol/(kg.bar)</p>
<ul>
<li>pO2 = 95 mmHg =&gt; dissolved O2 .. 0.14 mmol/kg at 37degC</li>
</ul>
</html>",
        revisions="<html>
<p><i>2015-2019</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
      Diagram(coordinateSystem(extent={{-160,-100},{160,100}})));
  end GasSolubility_NIST;

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
    Chemical.Components.Substance CO2_unbound_plasma(
      amountOfSubstance_start(displayUnit="mmol") = 0.001,
      substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
      use_mass_start=false) "Free dissolved CO2 in blood plasma"
      annotation (Placement(transformation(extent={{-90,-24},{-70,-4}})));
    Chemical.Components.GasSolubility O2_dissolutionP
      annotation (Placement(transformation(extent={{-34,44},{-14,64}})));

    Chemical.Sources.ExternalIdealGasSubstance O2_g_n1(
      substanceData=Chemical.Substances.Oxygen_gas(),
      PartialPressure=12665.626804425,
      TotalPressure=system.p_ambient)
      annotation (Placement(transformation(extent={{22,76},{42,96}})));
    Chemical.Components.Substance O2_unbound_plasma(
      amountOfSubstance_start(displayUnit="mmol") = 0.0001,
      substanceData=Chemical.Substances.Oxygen_aqueous(),
      use_mass_start=false) "Free dissolved O2 in blood plasma"
      annotation (Placement(transformation(extent={{-50,-26},{-30,-6}})));
    Chemical.Components.GasSolubility CO2_dissolutionE
      annotation (Placement(transformation(extent={{36,44},{56,64}})));

    Chemical.Sources.ExternalIdealGasSubstance CO2_g_n2(
      substanceData=Chemical.Substances.CarbonDioxide_gas(),
      TotalPressure=system.p_ambient,
      PartialPressure(displayUnit="mmHg") = 5332.8954966)
      annotation (Placement(transformation(extent={{-58,78},{-38,98}})));

    Chemical.Components.Substance CO2_unbound_erythrocyte(
      amountOfSubstance_start(displayUnit="mmol") = 0.001,
      substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
      use_mass_start=false) "Free dissolved CO2 in red cells"
      annotation (Placement(transformation(extent={{18,-32},{38,-12}})));

    Chemical.Components.GasSolubility O2_dissolutionE_NIST
      annotation (Placement(transformation(extent={{78,44},{98,64}})));
    Chemical.Components.Substance O2_unbound_erythrocyte_NIST(
      amountOfSubstance_start(displayUnit="mmol") = 0.0001,
      substanceData=Chemical.Substances.Oxygen_aqueous(),
      use_mass_start=false) "Free dissolved O2 in red cells"
      annotation (Placement(transformation(extent={{58,-32},{78,-12}})));
    Components.Substance water_plasma(substanceData=
          Chemical.Substances.Water_liquid(), mass_start=0.82)
      annotation (Placement(transformation(extent={{-40,-66},{-20,-46}})));
    Components.Substance water(substanceData=Chemical.Substances.Water_liquid(),
        mass_start=1)
      annotation (Placement(transformation(extent={{68,-68},{88,-48}})));
    inner Modelica.Fluid.System system(p_ambient(displayUnit="mmHg")=
        101325.0144354, T_ambient=310.15)
      annotation (Placement(transformation(extent={{-10,-96},{10,-76}})));
    Components.Substance water_plasma1(mass_start=0.18, substanceData=
          Chemical.Interfaces.Incompressible.SubstanceData(MolarWeight=1/0.627))
      annotation (Placement(transformation(extent={{-76,-66},{-56,-46}})));
  equation

  connect(CO2_g_n2.port_a, CO2_dissolutionP.gas_port) annotation (Line(
      points={{-38,88},{-26,88},{-26,72},{-68,72},{-68,64}},
      color={158,66,200},
      thickness=1));
  connect(CO2_g_n2.port_a, CO2_dissolutionE.gas_port) annotation (Line(
      points={{-38,88},{-26,88},{-26,72},{46,72},{46,64}},
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
    connect(water_plasma.solution, blood_plasma.solution) annotation (Line(points=
           {{-36,-66},{-36,-75.1},{-26.4,-75.1}}, color={127,127,0}));
    connect(water.solution, red_cells.solution) annotation (Line(points={{72,-68},
            {72,-77.08},{83.2,-77.08}}, color={127,127,0}));
    connect(water_plasma1.solution, blood_plasma.solution) annotation (Line(
          points={{-72,-66},{-72,-76},{-26.4,-75.1}}, color={127,127,0}));
    annotation (
      experiment(StopTime=1e-005),
      Documentation(info="<html>
<p>Demonstration of different blood gases solubility in erythrocytes and in plasma. The difference is governed by various amount of other substances in the solution. </p>
<p>Please note, that the total content of CO2 and O2 in blood plasma and erythrocytes must be determined by including bicarbonate and hemoglobin connected amounts. </p>
</html>",
        revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end GasSolubility;

  model GasSolubility_blood "Dissolution of gases in liquids"
     extends Modelica.Icons.Example;

    Chemical.Components.Solution blood_plasma
      annotation (Placement(transformation(extent={{-100,-76},{-8,14}})));
                                        //(amountOfSolution_start=52.3)
    //(amountOfSolution_start=39.7)
    Chemical.Components.GasSolubility CO2_dissolutionP
      annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
    Chemical.Components.Substance CO2_unbound_plasma(
      amountOfSubstance_start(displayUnit="mmol") = 0.001,
      substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
      use_mass_start=false) "Free dissolved CO2 in blood plasma"
      annotation (Placement(transformation(extent={{-90,-24},{-70,-4}})));
    Chemical.Components.GasSolubility O2_dissolutionP
      annotation (Placement(transformation(extent={{-34,44},{-14,64}})));

    Chemical.Sources.ExternalIdealGasSubstance O2_g_n1(
      substanceData=Chemical.Substances.Oxygen_gas(),
      PartialPressure=12665.626804425,
      TotalPressure=system.p_ambient)
      annotation (Placement(transformation(extent={{22,76},{42,96}})));
    Chemical.Components.Substance O2_unbound_plasma(
      amountOfSubstance_start(displayUnit="mmol") = 0.0001,
      substanceData=Chemical.Substances.Oxygen_aqueous(),
      use_mass_start=false) "Free dissolved O2 in blood plasma"
      annotation (Placement(transformation(extent={{-50,-26},{-30,-6}})));

    Chemical.Sources.ExternalIdealGasSubstance CO2_g_n2(
      substanceData=Chemical.Substances.CarbonDioxide_gas(),
      TotalPressure=system.p_ambient,
      PartialPressure(displayUnit="mmHg") = 5332.8954966)
      annotation (Placement(transformation(extent={{-58,78},{-38,98}})));

    Components.Substance water(substanceData=Chemical.Substances.Water_liquid(),
        mass_start=0.82)
      annotation (Placement(transformation(extent={{-40,-66},{-20,-46}})));
    inner Modelica.Fluid.System system(p_ambient(displayUnit="mmHg")=
        101325.0144354, T_ambient=310.15)
      annotation (Placement(transformation(extent={{-10,-96},{10,-76}})));
    Components.Substance others(mass_start=0.18, substanceData=
          Chemical.Interfaces.Incompressible.SubstanceData(MolarWeight=1/0.627))
      annotation (Placement(transformation(extent={{-76,-66},{-56,-46}})));
  equation

  connect(CO2_g_n2.port_a, CO2_dissolutionP.gas_port) annotation (Line(
      points={{-38,88},{-26,88},{-26,72},{-68,72},{-68,64}},
      color={158,66,200},
      thickness=1));
  connect(CO2_dissolutionP.liquid_port, CO2_unbound_plasma.port_a)
    annotation (Line(
      points={{-68,44},{-68,-14},{-70,-14}},
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
    connect(water.solution, blood_plasma.solution) annotation (Line(points={{-36,
            -66},{-36,-75.1},{-26.4,-75.1}}, color={127,127,0}));
    connect(others.solution, blood_plasma.solution) annotation (Line(points={{
            -72,-66},{-72,-76},{-26.4,-75.1}}, color={127,127,0}));
    annotation (
      experiment(StopTime=1e-005),
      Documentation(info="<html>
<p>Demonstration of different blood gases solubility in erythrocytes and in plasma. The difference is governed by various amount of other substances in the solution. </p>
<p>Please note, that the total content of CO2 and O2 in blood plasma and erythrocytes must be determined by including bicarbonate and hemoglobin connected amounts. </p>
</html>",
        revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end GasSolubility_blood;

  model EnzymeKinetics "Basic enzyme kinetics"
    extends Modelica.Icons.Example;

    Chemical.Components.Solution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
    Chemical.Components.Substance P(
      substanceData(DfG=-Modelica.Constants.R*298.15*50),
      use_mass_start=false,
      amountOfSubstance_start=1e-8)
      annotation (Placement(transformation(extent={{92,-12},{72,8}})));
    Chemical.Components.Substance S(use_mass_start=false, amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));

    parameter Modelica.Units.SI.AmountOfSubstance tE=1e-6
      "Total amount of enzyme";
       parameter Real k_cat(unit="1/s", displayUnit="1/min")= 1
      "Forward rate of second reaction";
    constant Modelica.Units.SI.Concentration Km=0.1
      "Michaelis constant = substrate concentration at rate of half Vmax";

    parameter Modelica.Units.SI.MolarFlowRate Vmax=1e-5*k_cat
      "Maximal molar flow";


    Chemical.Components.Substance ES(
      substanceData(DfG=-Modelica.Constants.R*298.15*log(2/Km)),
      use_mass_start=false,
      amountOfSubstance_start=tE/2)
      annotation (Placement(transformation(extent={{-12,-10},{8,10}})));
    Chemical.Components.Substance E(use_mass_start=false, amountOfSubstance_start=tE/2)
      annotation (Placement(transformation(extent={{-10,38},{10,58}})));
    Chemical.Components.Reaction chemicalReaction(
      nS=2,
      KC=Vmax/(2*Modelica.Constants.R*298.15*log(2)),
      nP=1)
      annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));

    Chemical.Components.Reaction chemicalReaction1(
      nP=2,
      KC=Vmax/(2*Modelica.Constants.R*298.15*(50 - log(2))),
      nS=1) annotation (Placement(transformation(extent={{24,-10},{44,10}})));

    Components.Substance liquidWater(substanceData=
          Chemical.Substances.Water_liquid(), mass_start=1)
      annotation (Placement(transformation(extent={{42,-80},{62,-60}})));
  equation
    //Michaelis-Menton: v=((E.q_out.conc + ES.q_out.conc)*k_cat)*S.concentration/(Km+S.concentration);
    connect(S.port_a, chemicalReaction.substrates[1]) annotation (Line(
        points={{-72,-4},{-56,-4},{-56,2},{-42,2}},
        color={158,66,200},
        thickness=1));
    connect(E.port_a, chemicalReaction.substrates[2]) annotation (Line(
        points={{10,48},{-52,48},{-52,-2},{-42,-2}},
        color={158,66,200},
        thickness=1));
    connect(E.port_a, chemicalReaction1.products[2]) annotation (Line(
        points={{10,48},{54,48},{54,-2},{44,-2}},
        color={158,66,200},
        thickness=1));
    connect(chemicalReaction1.products[1], P.port_a) annotation (Line(
        points={{44,2},{58,2},{58,-2},{72,-2}},
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
    connect(liquidWater.solution, solution.solution) annotation (Line(points={{
            46,-80},{46,-98},{60,-98}}, color={127,127,0}));
    connect(chemicalReaction.products[1], ES.port_a)
      annotation (Line(points={{-22,0},{8,0}}, color={158,66,200}));
    connect(ES.port_a, chemicalReaction1.substrates[1])
      annotation (Line(points={{8,0},{24,0}}, color={158,66,200}));
        annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
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
</html>"),
      experiment(StopTime=199000, __Dymola_Algorithm="Dassl"));
  end EnzymeKinetics;

  model ElectrochemicalCell
    "The electrochemical cell: Pt(s) | H2(g) | H+(aq), Cl-(aq) | AgCl(s) | Ag(s)"
   extends Modelica.Icons.Example;

    Chemical.Components.Solution cathode(ElectricGround=false)
      annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
    Chemical.Components.Substance Ag(
      substanceData=Chemical.Substances.Silver_solid(),
      use_mass_start=false,
      amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-72,-30},{-52,-10}})));
    Chemical.Components.Solution anode(ElectricGround=false)
      annotation (Placement(transformation(extent={{62,-50},{96,50}})));

    Chemical.Components.Solution solution1(ElectricGround=false)
      annotation (Placement(transformation(extent={{-30,-60},{38,6}})));

    Chemical.Components.Substance Cl(
      substanceData=Chemical.Substances.Chloride_aqueous(),
      use_mass_start=false,
      amountOfSubstance_start=12.39)
      annotation (Placement(transformation(extent={{0,-26},{-20,-6}})));
    Chemical.Components.Substance AgCl(
      substanceData=Chemical.Substances.SilverChloride_solid(),
      use_mass_start=false,
      amountOfSubstance_start=1e-8)
      annotation (Placement(transformation(extent={{-76,4},{-56,24}})));
    Chemical.Sources.ExternalIdealGasSubstance H2(
      substanceData=Chemical.Substances.Hydrogen_gas(),
      PartialPressure=100000,
      TotalPressure=100000)
      annotation (Placement(transformation(extent={{24,32},{44,52}})));
    Chemical.Components.Substance H(
      substanceData=Chemical.Substances.Proton_aqueous(),
      use_mass_start=false,
      amountOfSubstance_start=12.39)
      annotation (Placement(transformation(extent={{8,-26},{28,-6}})));
    Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
      annotation (Placement(transformation(extent={{-6,64},{14,84}})));
    Chemical.Components.Reaction electrodeReaction(
      nP=2,
      p={2,2},
      nS=1) annotation (Placement(transformation(
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
                               //(substanceData=Chemical.Examples.Substances.Electrone_solid())
  Chemical.Components.ElectronTransfer electrone1
      annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                //(substanceData=Chemical.Examples.Substances.Electrone_solid())
  Modelica.Electrical.Analog.Basic.Ground ground
    annotation (Placement(transformation(extent={{84,-84},{104,-64}})));
    Components.Substance liquidWater(substanceData=
          Chemical.Substances.Water_liquid(), mass_start=1)
      annotation (Placement(transformation(extent={{-6,-54},{14,-34}})));
  equation
    connect(Ag.port_a, electrodeReaction1.substrates[1]) annotation (Line(
        points={{-52,-20},{-42,-20},{-42,-10},{-38,-10}},
        color={158,66,200},
        thickness=1));
    connect(Cl.port_a, electrodeReaction1.substrates[2]) annotation (Line(
        points={{-20,-16},{-42,-16},{-42,-10}},
        color={158,66,200},
        thickness=1));
    connect(AgCl.port_a, electrodeReaction1.products[1]) annotation (Line(
        points={{-56,14},{-38,14},{-38,10}},
        color={158,66,200},
        thickness=1));
    connect(H.port_a, electrodeReaction.products[1]) annotation (Line(
        points={{28,-16},{50,-16},{50,-4}},
        color={158,66,200},
        thickness=1));
    connect(electrodeReaction.products[2], electrone1.port_a) annotation (Line(
        points={{54,-4},{54,-16},{68,-16}},
        color={158,66,200},
        thickness=1));
    connect(electrodeReaction1.products[2], electrone.port_a) annotation (Line(
        points={{-42,10},{-42,42},{-58,42}},
        color={158,66,200},
        thickness=1));
    connect(Cl.solution, solution1.solution) annotation (Line(
        points={{-4,-26},{-4,-30},{24.4,-30},{24.4,-59.34}},
        color={127,127,0}));
    connect(H.solution, solution1.solution) annotation (Line(points={{12,-26},{
            12,-30},{24.4,-30},{24.4,-59.34}},
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
    connect(liquidWater.solution, solution1.solution) annotation (Line(points={
            {-2,-54},{-2,-59.34},{24.4,-59.34}}, color={127,127,0}));
    connect(H2.port_a, electrodeReaction.substrates[1])
      annotation (Line(points={{44,42},{52,42},{52,16}}, color={158,66,200}));
    annotation (
    experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ElectrochemicalCell;

  model LeadAcidBattery
    "The electrochemical cell: PbSO4(s) | Pb(s) | HSO4-(aq) , H+(aq) | PbO2(s) | PbSO4(s) + 2 H2O"
   extends Modelica.Icons.Example;

    Chemical.Components.Solution anode(ElectricGround=false)
      annotation (Placement(transformation(extent={{24,-76},{58,32}})));

    Chemical.Components.Solution cathode(ElectricGround=false)
      annotation (Placement(transformation(extent={{-80,-78},{-46,30}})));

    Chemical.Components.Solution solution1(ElectricGround=false)
      annotation (Placement(transformation(extent={{-26,-80},{2,20}})));

    Chemical.Components.Substance Pb(
      substanceData=Chemical.Substances.Lead_solid(),
      use_mass_start=false,
      amountOfSubstance_start=50)
      annotation (Placement(transformation(extent={{50,-66},{30,-46}})));
    Chemical.Components.Substance HSO4(
      substanceData=Chemical.Substances.HydrogenSulfate_aqueous(),
      use_mass_start=false,
      amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-2,-70},{-22,-50}})));
    Chemical.Components.Substance PbSO4_(
      amountOfSubstance_start(displayUnit="mol") = 0.001,
      substanceData=Chemical.Substances.LeadSulfate_solid(),
      use_mass_start=false)
      annotation (Placement(transformation(extent={{52,-30},{32,-10}})));
    Chemical.Components.Substance H(
      substanceData=Chemical.Substances.Proton_aqueous(),
      use_mass_start=false,
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
    Chemical.Components.Substance PbO2(
      substanceData=Chemical.Substances.LeadDioxide_solid(),
      use_mass_start=false,
      amountOfSubstance_start=50) annotation (Placement(transformation(extent={{-10,-10},
              {10,10}}, origin={-60,-58})));
    Components.Substance H2O(mass_start(displayUnit="g") = 0.114, substanceData=
         Chemical.Substances.Water_liquid())
      annotation (Placement(transformation(extent={{-2,-8},{-22,12}})));
    Chemical.Components.Substance PbSO4(
      amountOfSubstance_start=0.001,
      substanceData=Chemical.Substances.LeadSulfate_solid(),
      use_mass_start=false) annotation (Placement(transformation(extent={{-10,-10},{10,
              10}}, origin={-60,6})));

  Modelica.Electrical.Analog.Basic.Ground ground
    annotation (Placement(transformation(extent={{16,30},{36,50}})));
  Modelica.Electrical.Analog.Basic.Resistor resistor(R=1)
    annotation (Placement(transformation(extent={{-14,40},{6,60}})));
  Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
    annotation (Placement(transformation(extent={{-56,40},{-36,60}})));

    Real density, molality, totalmolality, voltage;
    inner Modelica.Fluid.System system(T_ambient=299.15)
      annotation (Placement(transformation(extent={{62,64},{82,84}})));
  equation
    density = solution1.solution.m/solution1.solution.V;
    totalmolality = solution1.solution.n/((H2O.x*solution1.solution.n)*H2O.substanceData.MolarWeight);
    molality = HSO4.x*solution1.solution.n/((H2O.x*solution1.solution.n)*H2O.substanceData.MolarWeight);
    voltage = voltageSensor.v;

    connect(Pb.port_a, electrodeReaction1.substrates[1]) annotation (Line(
        points={{30,-56},{13.5,-56},{13.5,-26},{12,-26}},
        color={158,66,200},
        thickness=0.5));
    connect(HSO4.port_a, electrodeReaction1.substrates[2]) annotation (Line(
        points={{-22,-60},{16,-60},{16,-26}},
        color={158,66,200},
        thickness=0.5));
    connect(PbSO4_.port_a, electrodeReaction1.products[1]) annotation (Line(
        points={{32,-20},{26,-20},{26,-2},{16,-2},{16,-6},{11.3333,-6}},
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
        points={{-34,-4},{-34,6},{-50,6}},
        color={158,66,200},
        thickness=0.5));
    connect(electrodeReaction.products[2], H2O.port_a) annotation (Line(
        points={{-38,-4},{-38,2},{-22,2}},
        color={158,66,200},
        thickness=0.5));
    connect(PbO2.port_a, electrodeReaction.substrates[1]) annotation (Line(
        points={{-50,-58},{-36,-58},{-36,-24},{-33,-24}},
        color={158,66,200},
        thickness=0.5));
    connect(HSO4.port_a, electrodeReaction.substrates[2]) annotation (Line(
        points={{-22,-60},{-34,-60},{-34,-24},{-35,-24}},
        color={158,66,200},
        thickness=0.5));
    connect(H.port_a, electrodeReaction.substrates[3]) annotation (Line(
        points={{-22,-32},{-32,-32},{-32,-24},{-37,-24}},
        color={158,66,200},
        thickness=0.5));
    connect(electrone1.port_a, electrodeReaction.substrates[4]) annotation (Line(
        points={{-52,-28},{-38,-28},{-38,-24},{-39,-24}},
        color={158,66,200},
        thickness=0.5));
    connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(
        points={{-22,-32},{2,-32},{2,2},{12,2},{12,-6},{14,-6}},
        color={158,66,200},
        thickness=0.5));
    connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
        points={{30,12},{14,12},{14,-6},{16.6667,-6}},
        color={158,66,200},
        thickness=0.5));
  connect(Pb.solution, anode.solution) annotation (Line(
      points={{46,-66},{46,-74.92},{51.2,-74.92}},
      color={127,127,0}));
  connect(PbSO4_.solution, anode.solution) annotation (Line(
      points={{48,-30},{48,-74.92},{51.2,-74.92}},
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
    experiment(StopTime=47800, __Dymola_Algorithm="Dassl"),
                                Documentation(revisions=
                      "<html>
<p><i>2015-2018</i></p>
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
      Chemical.Components.Substance H3O(
        substanceData=Chemical.Substances.Hydronium_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{10,
                -10},{-10,10}}, origin={30,70})));
      Chemical.Components.Substance OH(
        substanceData=Chemical.Substances.Hydroxide_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{10,
                -10},{-10,10}}, origin={30,26})));
      Chemical.Components.Substance H2O(mass_start=1, substanceData=
            Chemical.Substances.Water_liquid()) annotation (Placement(
            transformation(extent={{-10,-10},{10,10}}, origin={-30,46})));
      Chemical.Components.Reaction waterDissociation(
        KC=1e-3,                                     nS=1, nP=2, s={2})
        annotation (Placement(transformation(extent={{-12,36},{8,56}})));
            Real pH, pH3O;
      Chemical.Components.Substance H_(
        substanceData=Chemical.Substances.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{10,
                -10},{-10,10}}, origin={28,-30})));
      Chemical.Components.Substance OH_(
        substanceData=Chemical.Substances.Hydroxide_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{10,
                -10},{-10,10}}, origin={28,-76})));
      Chemical.Components.Substance H2O_(mass_start=1, substanceData=
            Chemical.Substances.Water_liquid()) annotation (Placement(
            transformation(extent={{-10,-10},{10,10}}, origin={-32,-56})));
      Chemical.Components.Reaction waterDissociation_(
        KC=1e-3,                                      nS=1, nP=2)
        annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));

      inner Modelica.Fluid.System system(p_ambient=100000, T_ambient=298.15)
        annotation (Placement(transformation(extent={{-90,50},{-70,70}})));
    equation
      pH3O = -log10( H3O.a);

      pH = -log10( H_.a);

      connect(OH.port_a, waterDissociation.products[1]) annotation (Line(
          points={{20,26},{16,26},{16,48},{8,48}},
          color={158,66,200},
          thickness=1));
      connect(waterDissociation.products[2], H3O.port_a) annotation (Line(
          points={{8,44},{16,44},{16,70},{20,70}},
          color={158,66,200},
          thickness=1));
      connect(H2O.port_a, waterDissociation.substrates[1]) annotation (Line(
          points={{-20,46},{-12,46}},
          color={158,66,200},
          thickness=1));
      connect(OH_.port_a,waterDissociation_. products[1]) annotation (Line(
          points={{18,-76},{14,-76},{14,-54},{6,-54}},
          color={158,66,200},
          thickness=1));
      connect(waterDissociation_.products[2], H_.port_a) annotation (Line(
          points={{6,-58},{14,-58},{14,-30},{18,-30}},
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
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=1));
    end WaterSelfIonization;

    model CarbonDioxideInWater "CO2 as alone acid-base buffer"
      import Chemical;
        extends Modelica.Icons.Example;
        parameter Real KC = 1e-3;

      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,46}})));
      Chemical.Components.Substance HCO3(
        amountOfSubstance_start(displayUnit="mmol") = 1e-08,
        substanceData=Chemical.Substances.Bicarbonate_aqueous(),
        use_mass_start=false)
        annotation (Placement(transformation(extent={{-16,-4},{4,16}})));
      Chemical.Components.Reaction HendersonHasselbalch(
        KC=1e-4*KC,                                     nP=2, nS=2,
      useKineticsInput=false) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
        annotation (Placement(transformation(extent={{-48,-6},{-28,14}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
        substanceData=Chemical.Substances.CarbonDioxide_gas(),
        PartialPressure(displayUnit="mmHg") = 5332.8954966,
        TotalPressure(displayUnit="mmHg") = 101325.0144354) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-60,86})));
      Chemical.Components.Substance H(
        substanceData=Chemical.Substances.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-10,
                -10},{10,10}}, origin={10,-30})));
      Chemical.Components.GasSolubility gasSolubility(KC=KC)
        annotation (Placement(transformation(extent={{-70,36},{-50,56}})));
                                            /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
      Chemical.Components.Substance CO2_liquid(
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
        use_mass_start=false)
        annotation (Placement(transformation(extent={{-82,-6},{-62,14}})));
      Real pH;

      Chemical.Components.Substance liquidWater(substanceData=
            Chemical.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{-76,-50},{-56,-30}})));
      inner Modelica.Fluid.System system(T_ambient=310.15)
        annotation (Placement(transformation(extent={{48,64},{68,84}})));
      Chemical.Components.Substance OH(
        substanceData=Chemical.Substances.Hydroxide_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-10,
                -10},{10,10}}, origin={12,-72})));
      Chemical.Components.Reaction waterDissociation(
        KC=KC,
        nP=2,
        nS=1,
        useKineticsInput=false) "water dissociation" annotation (
          Placement(transformation(extent={{-44,-68},{-24,-48}})));
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
          points={{-28,6},{-22,6},{-22,-30},{20,-30}},
          color={158,66,200},
          thickness=1));
      connect(HendersonHasselbalch.products[2], HCO3.port_a) annotation (Line(
          points={{-28,2},{-12,2},{-12,6},{4,6}},
          color={158,66,200},
          thickness=1));
      connect(CO2_liquid.port_a, HendersonHasselbalch.substrates[2]) annotation (
          Line(
          points={{-62,4},{-62,2},{-48,2}},
          color={158,66,200},
          thickness=1));
      connect(CO2_liquid.solution, solution.solution) annotation (Line(
          points={{-78,-6},{-78,-98.54},{60,-98.54}},
          color={127,127,0}));
      connect(HCO3.solution, solution.solution) annotation (Line(points={{-12,-4},
            {-12,-98.54},{60,-98.54}},color={127,127,0}));
      connect(liquidWater.solution, solution.solution) annotation (Line(points={
              {-72,-50},{-72,-98.54},{60,-98.54}}, color={127,127,0}));
      connect(liquidWater.port_a, HendersonHasselbalch.substrates[1])
        annotation (Line(points={{-56,-40},{-54,-40},{-54,6},{-48,6}}, color={158,
              66,200}));
      connect(liquidWater.port_a, waterDissociation.substrates[1])
        annotation (Line(points={{-56,-40},{-50,-40},{-50,-58},{-44,-58}},
            color={158,66,200}));
      connect(waterDissociation.products[1], H.port_a) annotation (Line(
            points={{-24,-56},{-6,-56},{-6,-30},{20,-30}}, color={158,66,
              200}));
      connect(waterDissociation.products[2], OH.port_a) annotation (Line(
            points={{-24,-60},{-8,-60},{-8,-72},{22,-72}}, color={158,66,
              200}));
      connect(H.solution, solution.solution) annotation (Line(points={{4,
              -40},{-12,-40},{-12,-98.54},{60,-98.54}}, color={127,127,0}));
      connect(OH.solution, solution.solution) annotation (Line(points={{6,
              -82},{6,-98.54},{60,-98.54}}, color={127,127,0}));
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
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.02));
    end CarbonDioxideInWater;

    model Phosphate
      import Chemical;
        extends Modelica.Icons.Example;
        parameter Real KC = 1e-3;
      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-98,-100},{100,100}})));

      Chemical.Components.Substance H(
        substanceData=Chemical.Substances.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=10^(-7.4)) "hydrogen ions activity" annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={28,-14})));

      Chemical.Components.Substance H3PO4(
        substanceData=Chemical.Substances.PhosphoricAcid_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-08)
        annotation (Placement(transformation(extent={{-90,-58},{-70,-38}})));
      Chemical.Components.Substance H2PO4(
        substanceData=Chemical.Substances.DihydrogenPhosphate_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=0.0005)
        annotation (Placement(transformation(extent={{-40,-58},{-20,-38}})));
      Chemical.Components.Substance HPO4(
        substanceData=Chemical.Substances.HydrogenPhosphate_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=0.0006)
        annotation (Placement(transformation(extent={{16,-58},{36,-38}})));
      Chemical.Components.Substance PO4(
        substanceData=Chemical.Substances.Phosphate_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-08)
        annotation (Placement(transformation(extent={{92,-58},{72,-38}})));

      Chemical.Components.Reaction chemicalReaction(
        KC=KC,                                      nS=1, nP=2) "10^(-1.915 + 3)"
        annotation (Placement(transformation(extent={{-66,-58},{-46,-38}})));
      Chemical.Components.Reaction chemicalReaction1(
        KC=KC,                                       nS=1, nP=2) "10^(-6.66 + 3)"
        annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
      Chemical.Components.Reaction chemicalReaction2(
        KC=KC,                                       nS=1, nP=2) "10^(-11.78 + 3)"
        annotation (Placement(transformation(extent={{44,-58},{64,-38}})));

      Chemical.Components.Substance H2O(substanceData=
            Chemical.Substances.Water_liquid(), mass_start=1) annotation (
          Placement(transformation(extent={{10,-10},{-10,10}}, origin={86,14})));
      Real pH "acidity";
      Chemical.Components.Substance OH(
        substanceData=Chemical.Substances.Hydroxide_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=10^(-14 + 7.4)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}}, origin={28,32})));
      Chemical.Components.Reaction reaction(
        KC=KC,
        nS=2,
        nP=1) annotation (Placement(transformation(extent={{46,4},{66,24}})));
    equation
      pH = -log10( H.a);
      connect(H3PO4.port_a, chemicalReaction.substrates[1]) annotation (Line(
          points={{-70,-48},{-66,-48}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction.products[1], H2PO4.port_a) annotation (Line(
          points={{-46,-46},{-42,-46},{-42,-48},{-20,-48}},
          color={107,45,134},
          thickness=1));
      connect(H2PO4.port_a, chemicalReaction1.substrates[1]) annotation (Line(
          points={{-20,-48},{-14,-48}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction1.products[1], HPO4.port_a) annotation (Line(
          points={{6,-46},{16,-46},{16,-48},{36,-48}},
          color={107,45,134},
          thickness=1));
      connect(HPO4.port_a, chemicalReaction2.substrates[1]) annotation (Line(
          points={{36,-48},{44,-48}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction2.products[1], PO4.port_a) annotation (Line(
          points={{64,-46},{74,-46},{74,-48},{72,-48}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction.products[2], H.port_a) annotation (Line(
          points={{-46,-50},{-44,-50},{-44,-32},{38,-32},{38,-14}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction1.products[2], H.port_a) annotation (Line(
          points={{6,-50},{14,-50},{14,-32},{38,-32},{38,-14}},
          color={107,45,134},
          thickness=1));
      connect(chemicalReaction2.products[2], H.port_a) annotation (Line(
          points={{64,-50},{66,-50},{66,-32},{38,-32},{38,-14}},
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
        points={{92,4},{92,-98},{60.4,-98}},
        color={158,66,200}));
      connect(OH.solution, solution.solution) annotation (Line(points={{22,22},
              {22,-88},{60.4,-88},{60.4,-98}}, color={127,127,0}));
      connect(OH.port_a, reaction.substrates[1]) annotation (Line(points={{38,32},
              {42,32},{42,16},{46,16}}, color={158,66,200}));
      connect(H.port_a, reaction.substrates[2]) annotation (Line(points={{38,-14},
              {42,-14},{42,12},{46,12}}, color={158,66,200}));
      connect(reaction.products[1], H2O.port_a)
        annotation (Line(points={{66,14},{76,14}}, color={158,66,200}));
      annotation ( Documentation(info="<html>
</html>",    revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.05));
    end Phosphate;

    model CarbonDioxideInBlood
      import Chemical;
        extends Modelica.Icons.Example;

      parameter Real KC=1e-3;
      //e-6 "Slow down factor";
      Chemical.Components.Solution blood_erythrocytes(ElectricGround=false,
          temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-98},{100,-38}})));
      Chemical.Components.Solution blood_plasma(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,4},{100,56}})));

      Chemical.Components.Substance HCO3(
        substanceData=Chemical.Substances.Bicarbonate_blood(),
        use_mass_start=false,
        amountOfSubstance_start=0.024) annotation (Placement(transformation(extent={{
                10,-10},{-10,10}}, origin={18,24})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
        substanceData=Chemical.Substances.CarbonDioxide_gas(),
        TotalPressure(displayUnit="mmHg") = 101325.0144354,
        PartialPressure(displayUnit="mmHg") = 5332.8954966,
        usePartialPressureInput=true,
        Temperature=310.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-84,84})));
      Chemical.Components.GasSolubility gasSolubility(KC=KC)
        annotation (Placement(transformation(extent={{-94,48},{-74,68}})));

      Chemical.Components.Substance CO2(
        substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=0.00148) "Free dissolved CO2 in plasma"
        annotation (Placement(transformation(extent={{-88,28},{-68,48}})));
      Chemical.Components.Substance H2O(substanceData=
            Chemical.Substances.Water_liquid(), mass_start=51.6159/55.508)
        annotation (Placement(transformation(extent={{-60,12},{-40,32}})));
      Chemical.Components.Substance HCO3_E(
        substanceData=Chemical.Substances.Bicarbonate_blood(),
        use_mass_start=false,
        amountOfSubstance_start=0.0116)
        annotation (Placement(transformation(extent={{28,-60},{8,-40}})));
      Chemical.Components.Reaction HendersonHasselbalch1(nP=2, nS=2,
      KC=KC) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
        annotation (Placement(transformation(extent={{-26,-68},{-6,-48}})));
      Chemical.Components.Substance CO2_E(
        substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=0.0011) "Free dissolved CO2 in erythrocyte"
        annotation (Placement(transformation(extent={{-90,-82},{-70,-62}})));
      Chemical.Components.Substance H2O_E(substanceData=
            Chemical.Substances.Water_liquid(), mass_start=38.4008/55.508)
        annotation (Placement(transformation(extent={{-60,-62},{-40,-42}})));
      Chemical.Components.Substance Cl_E(
        substanceData=Chemical.Substances.Chloride_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=0.0499)
        annotation (Placement(transformation(extent={{68,-60},{48,-40}})));
      Chemical.Components.Substance Cl(
        substanceData=Chemical.Substances.Chloride_aqueous(),
        use_mass_start=false,
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
        substanceData=Chemical.Substances.Proton_aqueous(),
        BufferValue=0.063,
        a_start=10^(-7.2))
        annotation (Placement(transformation(extent={{48,-84},{30,-66}})));
      Modelica.Blocks.Sources.ContinuousClock clock(offset=5000)
        annotation (Placement(transformation(extent={{-54,62},{-34,82}})));
      Chemical.Components.Substance others_E(
        substanceData=Chemical.Interfaces.Incompressible.SubstanceData(
                density=(1.045 - 0.695523)*1000/(1 - 0.697583),
                References={"erythrocyte intracellular fluid density 1045kg/m3"},
                MolarWeight=(1.045 - 0.695523)/(38.7*(1 - 0.994648) - 0.0499 - 0.0116
               - 0.00123)),
        use_mass_start=false,
        amountOfSubstance_start=0.1444)
        annotation (Placement(transformation(extent={{68,-88},{88,-68}})));
      Chemical.Components.Substance others_P(
        substanceData=Chemical.Interfaces.Incompressible.SubstanceData(
                References={
              "to reach plasma density 1024 kg/m3 and plasma volume 1 liter"},
                density=(1.024 - 0.933373)*1000/(1 - 0.936137),
                MolarWeight=(1.024 - 0.933373)/(51.8*(1 - 0.994648) - 0.103 - 0.024 -
              0.0017)),
        use_mass_start=false,
        amountOfSubstance_start=0.1487)
        annotation (Placement(transformation(extent={{70,14},{90,34}})));
      Chemical.Components.Diffusion diffusion annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-66,-16})));
      Chemical.Sources.Buffer H(
        substanceData=Chemical.Substances.Proton_aqueous(),
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
          points={{-6,-56},{2,-56},{2,-50},{8,-50}},
          color={107,45,134},
          thickness=0.5));
    connect(CO2_E.port_a, HendersonHasselbalch1.substrates[1]) annotation (
        Line(
        points={{-70,-72},{-36,-72},{-36,-56},{-26,-56}},
        color={107,45,134},
        thickness=0.5));
      connect(H2O_E.port_a, HendersonHasselbalch1.substrates[2]) annotation (Line(
          points={{-40,-52},{-34,-52},{-34,-60},{-26,-60}},
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
        points={{30,-75},{4,-75},{4,-60},{-6,-60}},
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
          points={{-68,38},{-48,38},{-48,34},{-26,34}},
          color={158,66,200},
          thickness=0.5));
      connect(H2O.port_a, HendersonHasselbalch2.substrates[1]) annotation (Line(
          points={{-40,22},{-34,22},{-34,38},{-26,38}},
          color={158,66,200},
          thickness=0.5));
      connect(HendersonHasselbalch2.products[1], HCO3.port_a) annotation (Line(
          points={{-6,38},{2,38},{2,24},{8,24}},
          color={158,66,200},
          thickness=0.5));
      connect(HendersonHasselbalch2.products[2], H.port_a) annotation (Line(
          points={{-6,34},{2,34},{2,48},{20,48},{20,47}},
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
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(
        StopTime=3),
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
            Chemical.Substances.Proton_aqueous(), MoleFraction=10^(-7.1))
        annotation (Placement(transformation(extent={{0,-46},{20,-26}})));
      Components.Substance liquidWater(substanceData=
            Chemical.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{40,-80},{60,-60}})));
    equation
      connect(buffer.solution, simpleSolution.solution) annotation (Line(
          points={{-46,4},{-26,4},{-26,-98},{56,-98}},
          color={127,127,0}));
      connect(externalMoleFraction.port_a, buffer.port_a) annotation (Line(
          points={{20,-36},{40,-36},{40,10},{-30,10},{-30,14}},
          color={158,66,200}));
      connect(liquidWater.solution, simpleSolution.solution)
        annotation (Line(points={{44,-80},{44,-98},{56,-98}}, color={127,127,0}));
      annotation (                experiment(StopTime=0.05));
    end AcidBaseBufferTest;

    package Dev
      model RedCellMembrane
       // import Chemical;
        extends Modelica.Icons.Example;

        parameter Real KC=1;
        //e-6 "Slow down factor";
        Chemical.Components.Solution blood_erythrocytes(ElectricGround=false)
          annotation (Placement(transformation(extent={{-180,-100},{180,-10}})));
        Chemical.Components.Solution blood_plasma annotation (Placement(
              transformation(extent={{-180,12},{180,100}})));

        Chemical.Components.Substance HCO3(
          substanceData=Chemical.Substances.Bicarbonate_blood(),
          use_mass_start=false,
          amountOfSubstance_start=0.024) annotation (Placement(transformation(extent={
                  {-10,-10},{10,10}}, origin={-18,30})));

        Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
            mass_start=51.8*0.994648/55.508)
          annotation (Placement(transformation(extent={{-146,44},{-166,64}})));
        Chemical.Components.Substance HCO3_E(
          substanceData=Chemical.Substances.Bicarbonate_blood(),
          use_mass_start=false,
          amountOfSubstance_start=0.0116)
          annotation (Placement(transformation(extent={{-28,-38},{-8,-18}})));
        Components.Substance H2O_E(substanceData=
              Chemical.Substances.Water_liquid(), mass_start=38.7*0.994648/
              55.508) annotation (Placement(transformation(extent={{-144,-38},{
                  -164,-18}})));
        Chemical.Components.Substance Cl_E(
          substanceData=Chemical.Substances.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.0499)
          annotation (Placement(transformation(extent={{-4,-38},{16,-18}})));
        Chemical.Components.Substance Cl(
          substanceData=Chemical.Substances.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.103)
          annotation (Placement(transformation(extent={{-4,20},{16,40}})));
        Chemical.Components.Substance albumin(
          substanceData(
            MolarWeight=66.463,
            z=-17,
            density=1080),
          use_mass_start=false,
          amountOfSubstance_start=0.0007)
          annotation (Placement(transformation(extent={{112,76},{92,96}})));
        Real pH_e,pH_p;

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
        Chemical.Components.Substance permeableUncharged(use_mass_start=false,
            amountOfSubstance_start=0.0118)
          annotation (Placement(transformation(extent={{166,20},{146,40}})));
        Chemical.Components.Substance permeableUncharged_E(
          substanceData(MolarWeight=0.1),
          use_mass_start=false,
          amountOfSubstance_start=0.00903)
          annotation (Placement(transformation(extent={{164,-38},{144,-18}})));
        Chemical.Components.Substance chargedImpermeable_E(
          substanceData(MolarWeight=1),
          use_mass_start=false,
          amountOfSubstance_start=0.0165)
          annotation (Placement(transformation(extent={{144,-62},{164,-42}})));
        Chemical.Components.Membrane leak(useKineticsInput=false, KC=KC)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={140,0})));
        Chemical.Components.Substance Lac_E(
          substanceData=Chemical.Substances.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00062)
          annotation (Placement(transformation(extent={{56,-38},{76,-18}})));
        Chemical.Components.Substance Lac(
          substanceData=Chemical.Substances.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00131)
          annotation (Placement(transformation(extent={{56,20},{76,40}})));
        Chemical.Components.Membrane MCT_(useKineticsInput=false, KC=KC)
          "Monocarboxylate transporters" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={78,0})));
        Chemical.Components.Substance H_E(
          substanceData=Chemical.Substances.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=10^(-7.2)) "H+"
          annotation (Placement(transformation(extent={{30,-38},{50,-18}})));
        Chemical.Components.Substance H(
          substanceData=Chemical.Substances.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=10^(-7.4)) "H+ in plasma"
          annotation (Placement(transformation(extent={{30,20},{50,40}})));
        Chemical.Components.Membrane MCT(useKineticsInput=false, KC=KC)
          "Monocarboxylate transporters" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={52,0})));
        Chemical.Components.Substance CO2(
          substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00167) "free dissolved unbound CO2"
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        Chemical.Components.Substance CO2_E(
          substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00125) "free dissolved unbound CO2"
          annotation (Placement(transformation(extent={{-58,-38},{-38,-18}})));
        Chemical.Components.Membrane freeCO2(KC=KC) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-38,2})));
        Chemical.Components.Substance O2(
          substanceData=Chemical.Substances.Oxygen_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.000167) "free dissolved undound oxygen"
          annotation (Placement(transformation(extent={{96,20},{116,40}})));
        Chemical.Components.Membrane freeO2(KC=KC) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={118,0})));
        Chemical.Components.Substance O2_E(
          substanceData=Chemical.Substances.Oxygen_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.000125) "free dissolved undound O2"
          annotation (Placement(transformation(extent={{96,-38},{116,-18}})));
        Chemical.Components.Substance K(
          substanceData=Chemical.Substances.Potassium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.004)
          annotation (Placement(transformation(extent={{-100,20},{-120,40}})));
        Chemical.Components.Substance Na(
          substanceData=Chemical.Substances.Sodium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.138)
          annotation (Placement(transformation(extent={{-124,20},{-144,40}})));
        Chemical.Components.Substance Na_E(
          substanceData=Chemical.Substances.Sodium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.007)
          annotation (Placement(transformation(extent={{-118,-38},{-138,-18}})));
        Chemical.Components.Substance K_E(
          substanceData=Chemical.Substances.Potassium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.096)
          annotation (Placement(transformation(extent={{-112,-38},{-92,-18}})));
        Chemical.Components.Substance H2PO4_E(
          substanceData=Chemical.Substances.DihydrogenPhosphate_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.000175)
          annotation (Placement(transformation(extent={{-84,-38},{-64,-18}})));
        Chemical.Components.Substance ADP_E(
          substanceData(z=-3),
          use_mass_start=false,
          amountOfSubstance_start=9.6e-05)
          annotation (Placement(transformation(extent={{-114,-62},{-94,-42}})));
        Chemical.Components.Substance ATP_E(
          substanceData(
            z=-4,
            DfH=16700,
            DfG=30500,
            References={
                "http://www.wiley.com/college/pratt/0471393878/student/review/thermodynamics/7_relationship.html"}),
          use_mass_start=false,
          amountOfSubstance_start=0.00128)
          annotation (Placement(transformation(extent={{-146,-62},{-166,-42}})));

        Chemical.Components.Substance HPO4_E(
          substanceData=Chemical.Substances.HydrogenPhosphate_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.000495)
          annotation (Placement(transformation(extent={{-84,-62},{-64,-42}})));
        Chemical.Components.Substance globulins(
          substanceData(
            MolarWeight=34,
            z=-2.43,
            density=1080),
          use_mass_start=false,
          amountOfSubstance_start=0.00082)
          annotation (Placement(transformation(extent={{150,76},{130,96}})));
        Chemical.Components.Substance Ca(
          substanceData=Chemical.Substances.Calcium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00175) "Ca2+"
          annotation (Placement(transformation(extent={{-78,20},{-98,40}})));
        Chemical.Components.Substance Mg(
          substanceData=Chemical.Substances.Magnesium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00108) "Mg2+"
          annotation (Placement(transformation(extent={{-112,-84},{-92,-64}})));
        Chemical.Components.Substance DPG(
          substanceData(
            MolarWeight=0.266,
            z=-2.2,
            density=1000),
          use_mass_start=false,
          amountOfSubstance_start=0.0051)
          annotation (Placement(transformation(extent={{128,-94},{108,-74}})));
        Chemical.Components.Substance GSH(
          substanceData(
            MolarWeight=0.2,
            z=-1,
            density=1000),
          use_mass_start=false,
          amountOfSubstance_start=0.00223)
          annotation (Placement(transformation(extent={{164,-94},{144,-74}})));
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
        pH_p = -log10(H.a);
        pH_e = -log10(H_E.a);
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
        connect(HPO4_E.solution, H2PO4_E.solution) annotation (Line(
            points={{-80,-62},{-110,-62},{-110,-38},{-80,-38}},
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
        connect(CO2.port_a, HendersonHasselbalch.substrates[2]) annotation (Line(
            points={{-40,30},{-38,30},{-38,56},{-34,56}},
            color={158,66,200},
            thickness=1));
        connect(HCO3.port_a, HendersonHasselbalch.products[2]) annotation (Line(
            points={{-8,30},{-8,30},{-8,56},{-14,56}},
            color={158,66,200},
            thickness=1));
        connect(HendersonHasselbalch.substrates[1], H2O.port_a) annotation (Line(
            points={{-34,52},{-166,52},{-166,54}},
            color={158,66,200},
            thickness=1));
        connect(HendersonHasselbalch.products[1], H.port_a) annotation (Line(
            points={{-14,52},{18,52},{50,52},{50,30}},
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
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=1e-008),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},
                  {180,100}}), graphics));
      end RedCellMembrane;

      model NaKATPase
        Components.Solution ICF(ElectricGround=false)
          annotation (Placement(transformation(extent={{-100,-100},{-20,98}})));
        Components.Solution ECF
          annotation (Placement(transformation(extent={{28,-100},{100,100}})));
        Components.Substance water_ICF(substanceData=
              Chemical.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-82,70},{-62,90}})));
        Components.Substance water_EFC(substanceData=
              Chemical.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{58,56},{38,76}})));
        Components.Reaction NaKATPase(
          nS=4,
          nP=4,
          s={3,2,1,1},
          p={3,2,1,1}) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={0,8})));
        Components.Substance Na_ICF(
          amountOfSubstance_start(displayUnit="mmol") = 0.01,
          substanceData=Chemical.Substances.Sodium_aqueous(),
          use_mass_start=false)
          annotation (Placement(transformation(extent={{-82,38},{-62,58}})));
        Components.Substance K_ICF(
          amountOfSubstance_start(displayUnit="mmol") = 0.14,
          substanceData=Chemical.Substances.Potassium_aqueous(),
          use_mass_start=false)
          annotation (Placement(transformation(extent={{-82,10},{-62,30}})));
        Components.Substance ATP4_ICF(
          amountOfSubstance_start(displayUnit="mmol") = 0.0038,
          substanceData=Chemical.Substances.ATP4_aqueous(),
          use_mass_start=false)
          annotation (Placement(transformation(extent={{-82,-16},{-62,4}})));
        Components.Substance Na_ECF(
          amountOfSubstance_start(displayUnit="mmol") = 0.14,
          substanceData=Chemical.Substances.Sodium_aqueous(),
          use_mass_start=false)
          annotation (Placement(transformation(extent={{58,-26},{38,-6}})));
        Components.Substance K_ECF(
          amountOfSubstance_start(displayUnit="mmol") = 0.005,
          substanceData=Chemical.Substances.Potassium_aqueous(),
          use_mass_start=false)
          annotation (Placement(transformation(extent={{58,28},{38,48}})));
        Components.Substance ADP3(
          substanceData=Chemical.Substances.ADP3_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.005/150)
          annotation (Placement(transformation(extent={{-82,-42},{-62,-22}})));
        Components.Substance H2PO4_ICF(
          amountOfSubstance_start(displayUnit="mmol") = 0.036,
          substanceData=Chemical.Substances.DihydrogenPhosphate_aqueous(),
          use_mass_start=false)
          annotation (Placement(transformation(extent={{-82,-68},{-62,-48}})));
        inner Modelica.Fluid.System system(T_ambient=310.15)
          annotation (Placement(transformation(extent={{-6,-90},{14,-70}})));
        Components.Substance Cl_ICF(
          amountOfSubstance_start(displayUnit="mmol") = 0.0987,
          substanceData=Chemical.Substances.Chloride_aqueous(),
          use_mass_start=false)
          annotation (Placement(transformation(extent={{-82,-92},{-62,-72}})));
        Components.Substance Cl_ECF(
          amountOfSubstance_start(displayUnit="mmol") = 0.145,
          substanceData=Chemical.Substances.Chloride_aqueous(),
          use_mass_start=false)
          annotation (Placement(transformation(extent={{50,-72},{70,-52}})));
      equation
        connect(water_ICF.solution, ICF.solution) annotation (Line(points={{-78,70},
                {-78,-100},{-36,-100},{-36,-98.02}}, color={127,127,0}));
        connect(water_EFC.solution, ECF.solution) annotation (Line(points={{54,56},
                {54,-98},{85.6,-98}}, color={127,127,0}));
        connect(water_EFC.solution, Na_ECF.solution)
          annotation (Line(points={{54,56},{54,-26}}, color={127,127,0}));
        connect(water_EFC.solution, K_ECF.solution)
          annotation (Line(points={{54,56},{54,28}}, color={127,127,0}));
        connect(water_ICF.solution, Na_ICF.solution)
          annotation (Line(points={{-78,70},{-78,38}}, color={127,127,0}));
        connect(water_ICF.solution, K_ICF.solution)
          annotation (Line(points={{-78,70},{-78,10}}, color={127,127,0}));
        connect(water_ICF.solution, ATP4_ICF.solution)
          annotation (Line(points={{-78,70},{-78,-16}}, color={127,127,0}));
        connect(water_ICF.solution, ADP3.solution)
          annotation (Line(points={{-78,70},{-78,-42}}, color={127,127,0}));
        connect(water_ICF.solution, H2PO4_ICF.solution)
          annotation (Line(points={{-78,70},{-78,-68}}, color={127,127,0}));
        connect(Na_ICF.port_a, NaKATPase.substrates[1])
          annotation (Line(points={{-62,48},{-3,48},{-3,18}}, color={158,66,200}));
        connect(K_ECF.port_a, NaKATPase.substrates[2])
          annotation (Line(points={{38,38},{-1,38},{-1,18}}, color={158,66,200}));
        connect(NaKATPase.products[1], Na_ECF.port_a) annotation (Line(points={{-3,
                -2},{-2,-2},{-2,-16},{38,-16}}, color={158,66,200}));
        connect(NaKATPase.products[2], K_ICF.port_a) annotation (Line(points={{-1,
                -2},{-1,-26},{-12,-26},{-12,20},{-62,20}}, color={158,66,200}));
        connect(ATP4_ICF.port_a, NaKATPase.substrates[3]) annotation (Line(points={
                {-62,-6},{-36,-6},{-36,32},{1,32},{1,18}}, color={158,66,200}));
        connect(water_ICF.port_a, NaKATPase.substrates[4])
          annotation (Line(points={{-62,80},{3,80},{3,18}}, color={158,66,200}));
        connect(NaKATPase.products[3], ADP3.port_a) annotation (Line(points={{1,-2},
                {2,-2},{2,-32},{-62,-32}}, color={158,66,200}));
        connect(H2PO4_ICF.port_a, NaKATPase.products[4])
          annotation (Line(points={{-62,-58},{3,-58},{3,-2}}, color={158,66,200}));
        connect(water_ICF.solution, Cl_ICF.solution)
          annotation (Line(points={{-78,70},{-78,-92}}, color={127,127,0}));
        connect(water_EFC.solution, Cl_ECF.solution)
          annotation (Line(points={{54,56},{54,-72}}, color={127,127,0}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end NaKATPase;
    end Dev;

    model AlbuminTitration "Figge-Fencl model (22. Dec. 2007)"
      extends Modelica.Icons.Example;

      Chemical.Components.Solution solution(redeclare package stateOfMatter =
            Interfaces.Incompressible)
        annotation (Placement(transformation(extent={{-104,-100},{96,100}})));

      constant Integer n=218 "Number of weak acid group in albumin molecule";
      constant Real pKAs[n]=cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
        "acid dissociation constants";
      constant Real K[n]=fill(10.0, n) .^ (-pKAs);
      constant Real DfG[n]= Modelica.Constants.R*(298.15)*log(K);

      Chemical.Components.Substance A[n](
        substanceData(each z=-1),
        each use_mass_start=false,
        each amountOfSubstance_start=0.00033) "deprotonated acid groups"
        annotation (Placement(transformation(extent={{26,-16},{6,4}})));
      Chemical.Components.Reaction react[n](
        each KC=1e-9,                            each nS=1, each nP=2)
        annotation (Placement(transformation(extent={{-44,-2},{-24,18}})));

      Chemical.Components.Substance HA[n](
        substanceData(DfG=DfG),
        each use_mass_start=false,
        each amountOfSubstance_start=0.00033) "protonated acid groups"
        annotation (Placement(transformation(extent={{-78,-2},{-58,18}})));

      Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
          mass_start=1) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={62,-68})));
      Sources.ExternalMoleFraction H(substanceData=
            Chemical.Substances.Proton_aqueous(), MoleFraction=10^(-7.4))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={18,42})));
    equation
      connect(react.products[1], A.port_a) annotation (Line(
          points={{-24,10},{-12,10},{-12,-6},{6,-6}},
          color={107,45,134},
          thickness=1));
      for i in 1:n loop
        connect(HA[i].solution, solution.solution) annotation (Line(
          points={{-74,-2},{-74,-86},{56,-86},{56,-98}},
          color={127,127,0}));
        connect(A[i].solution, solution.solution) annotation (Line(
          points={{22,-16},{22,-86},{56,-86},{56,-98}},
          color={127,127,0}));
        connect(H.port_a, react[i].products[2]) annotation (Line(points={{8,42},{-8,42},{
              -8,6},{-24,6}}, color={158,66,200}));
      end for;
      connect(HA.port_a, react.substrates[1]) annotation (Line(
          points={{-58,8},{-44,8}},
          color={107,45,134},
          thickness=1));

      connect(solution.solution, H2O.solution) annotation (Line(
        points={{56,-98},{56,-78}},
        color={127,127,0}));


      annotation ( Documentation(revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>The titration slope der(pH)/der(SID)=185 1/(mol/L) at pH=7.4 and tAlb=0.66 mmol/l.</p>
<p>Data and model is described in</p>
<p><font style=\"color: #222222; \">Jame Figge: Role of non-volatile weak acids (albumin, phosphate and citrate). In: Stewart&apos;s Textbook of Acid-Base, 2nd Edition, John A. Kellum, Paul WG Elbers editors, &nbsp;AcidBase org, 2009, pp. 216-232.</font></p>
</html>"));
    end AlbuminTitration;

  end AcidBase;

  package Hemoglobin "Hemoglobin blood gases binding"
    model Allosteric_Hemoglobin_MWC "Monod,Wyman,Changeux (1965)"
      extends Modelica.Icons.Example;

      constant Modelica.Units.SI.AmountOfSubstance THb=0.001
        "Total amount of hemoglobin";

      constant Modelica.Units.SI.Temperature T=298.15 "Base Temperature";
      constant Real RT=Modelica.Constants.R*T;


      constant Modelica.Units.SI.Volume OneLiter=0.001;

      constant Real L=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      constant Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      constant Modelica.Units.SI.Concentration KR=0.000671946*(55.508/
          38.7)
        "Oxygen dissociation coefficient on relaxed(R) hemoglobin subunit";

      constant Real KRx=(KR*OneLiter)
        "Mole fraction based KR";

      //Relative Gibbs formation energies of the substances in the system:
      constant Modelica.Units.SI.MolarEnergy GO2aq=-RT*log(0.0013);
      constant Modelica.Units.SI.MolarEnergy GR0=0;
      constant Modelica.Units.SI.MolarEnergy GT0=GR0 - RT*log(L);
      constant Modelica.Units.SI.MolarEnergy GR1=GR0 + GO2aq + RT*log(KRx
          /4);
      constant Modelica.Units.SI.MolarEnergy GT1=GR1 - RT*log(c*L);
      constant Modelica.Units.SI.MolarEnergy GR2=GR1 + GO2aq + RT*log(KRx
          /(3/2));
      constant Modelica.Units.SI.MolarEnergy GT2=GR2 - RT*log(c^2*L);
      constant Modelica.Units.SI.MolarEnergy GR3=GR2 + GO2aq + RT*log(KRx
          /(2/3));
      constant Modelica.Units.SI.MolarEnergy GT3=GR3 - RT*log(c^3*L);
      constant Modelica.Units.SI.MolarEnergy GR4=GR3 + GO2aq + RT*log(KRx
          *4);
      constant Modelica.Units.SI.MolarEnergy GT4=GR4 - RT*log(c^4*L);
                                      //*0.018),
      parameter Real KC = 0.0001 "Slow down factor";

      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-72,-102},{94,124}})));

      Chemical.Components.Substance oxygen_unbound(
        substanceData(DfG=GO2aq),
        use_mass_start=false,
        amountOfSubstance_start=1e-7)
        annotation (Placement(transformation(extent={{-62,-46},{-42,-26}})));

      Chemical.Components.Substance T0(
        substanceData(DfG=GT0),
        use_mass_start=false,
        amountOfSubstance_start=(THb))
        annotation (Placement(transformation(extent={{34,78},{54,98}})));

      Chemical.Components.Substance T1(
        substanceData(DfG=GT1),
        use_mass_start=false,
        amountOfSubstance_start=(THb*1e-4))
        annotation (Placement(transformation(extent={{34,36},{54,56}})));

      Chemical.Components.Substance T2(
        substanceData(DfG=GT2),
        use_mass_start=false,
        amountOfSubstance_start=THb*1e-8)
        annotation (Placement(transformation(extent={{34,-10},{54,10}})));

      Chemical.Components.Substance R1(
        substanceData(DfG=GR1),
        use_mass_start=false,
        amountOfSubstance_start=THb*1e-8)
        annotation (Placement(transformation(extent={{-20,36},{0,56}})));

      Chemical.Components.Substance R2(
        substanceData(DfG=GR2),
        use_mass_start=false,
        amountOfSubstance_start=THb*1e-10)
        annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

      Chemical.Components.Substance T3(
        substanceData(DfG=GT3),
        use_mass_start=false,
        amountOfSubstance_start=THb*1e-12)
        annotation (Placement(transformation(extent={{34,-54},{54,-34}})));

      Chemical.Components.Substance R3(
        substanceData(DfG=GR3),
        use_mass_start=false,
        amountOfSubstance_start=THb*1e-12)
        annotation (Placement(transformation(extent={{-20,-54},{0,-34}})));

      Chemical.Components.Substance T4(
        substanceData(DfG=GT4),
        use_mass_start=false,
        amountOfSubstance_start=THb*1e-17)
        annotation (Placement(transformation(extent={{34,-92},{54,-72}})));

      Chemical.Components.Substance R4(
        substanceData(DfG=GR4),
        use_mass_start=false,
        amountOfSubstance_start=THb*1e-14)
        annotation (Placement(transformation(extent={{-20,-92},{0,-72}})));

      Chemical.Components.Substance R0(
        substanceData(DfG=GR0),
        use_mass_start=false,
        amountOfSubstance_start=THb*1e-7)
        annotation (Placement(transformation(extent={{-20,78},{0,98}})));

      Chemical.Components.Reaction quaternaryForm(
        nS=1,
        nP=1,
        KC=KC) annotation (Placement(transformation(extent={{4,78},{24,98}})));
      Chemical.Components.Reaction oxyR1(
        nS=1,
        nP=2,
        KC=KC) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-8,64})));
      Chemical.Components.Reaction oxyT1(
        nS=1,
        nP=2,
        KC=KC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={44,64})));
      Chemical.Components.Reaction oxyR2(
        nS=1,
        nP=2,
        KC=KC) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-10,22})));
      Chemical.Components.Reaction oxyR3(
        nS=1,
        nP=2,
        KC=KC) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-10,-24})));
      Chemical.Components.Reaction oxyR4(
        nS=1,
        nP=2,
        KC=KC) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-10,-66})));
      Chemical.Components.Reaction oxyT2(
        nS=1,
        nP=2,
        KC=KC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={44,22})));
      Chemical.Components.Reaction oxyT3(
        nS=1,
        nP=2,
        KC=KC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={44,-24})));
      Chemical.Components.Reaction oxyT4(
        nS=1,
        nP=2,
        KC=KC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={44,-66})));
      Chemical.Components.Reaction quaternaryForm1(
        nS=1,
        nP=1,
        KC=KC) annotation (Placement(transformation(extent={{8,36},{28,56}})));
      Chemical.Components.Reaction quaternaryForm2(
        nS=1,
        nP=1,
        KC=KC)
        annotation (Placement(transformation(extent={{8,-10},{28,10}})));
      Chemical.Components.Reaction quaternaryForm3(
        nS=1,
        nP=1,
        KC=KC)
        annotation (Placement(transformation(extent={{8,-54},{28,-34}})));
      Chemical.Components.Reaction quaternaryForm4(
        nS=1,
        nP=1,
        KC=KC)
        annotation (Placement(transformation(extent={{10,-92},{30,-72}})));

      Modelica.Blocks.Sources.ContinuousClock clock(offset=10)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-84,62})));
      Chemical.Sources.ExternalIdealGasSubstance O2_in_air(
        usePartialPressureInput=true,
        TotalPressure(displayUnit="kPa") = 101325.0144354,
        substanceData=Chemical.Substances.Oxygen_gas(),
        PartialPressure(displayUnit="kPa") = 3733) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-84,22})));

      Chemical.Components.GasSolubility gasSolubility(KC=KC)
        annotation (Placement(transformation(extent={{-94,-16},{-74,4}})));

      Real sO2;
      Components.Substance substance(substanceData=Substances.Water_liquid(),
                                                mass_start=1)
        annotation (Placement(transformation(extent={{66,-72},{86,-52}})));
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
          points={{0,46},{0,32},{-8,32}},
          color={107,45,134}));
      connect(oxyR2.substrates[1],R2. port_a) annotation (Line(
          points={{-10,12},{-10,0},{0,0}},
          color={107,45,134}));
      connect(oxyR3.substrates[1],R3. port_a) annotation (Line(
          points={{-10,-34},{-10,-44},{0,-44}},
          color={107,45,134}));
      connect(oxyR3.products[1],R2. port_a) annotation (Line(
          points={{-8,-14},{-8,-7},{0,-7},{0,0}},
          color={107,45,134}));
      connect(R3.port_a,oxyR4. products[1]) annotation (Line(
          points={{0,-44},{0,-56},{-8,-56}},
          color={107,45,134}));
      connect(oxyR4.substrates[1],R4. port_a) annotation (Line(
          points={{-10,-76},{-10,-82},{0,-82}},
          color={107,45,134}));
      connect(oxyT1.products[1],T0. port_a) annotation (Line(
          points={{42,74},{42,88},{54,88}},
          color={107,45,134}));
      connect(oxyT1.substrates[1],T1. port_a) annotation (Line(
          points={{44,54},{44,46},{54,46}},
          color={107,45,134}));
      connect(T1.port_a,oxyT2. products[1]) annotation (Line(
          points={{54,46},{54,32},{42,32}},
          color={107,45,134}));
      connect(oxyT3.substrates[1],T3. port_a) annotation (Line(
          points={{44,-34},{44,-44},{54,-44}},
          color={107,45,134}));
      connect(T3.port_a,oxyT4. products[1]) annotation (Line(
          points={{54,-44},{54,-56},{42,-56}},
          color={107,45,134}));
      connect(oxyT4.substrates[1],T4. port_a) annotation (Line(
          points={{44,-76},{44,-82},{54,-82}},
          color={107,45,134}));
      connect(R0.port_a,quaternaryForm. substrates[1]) annotation (Line(
          points={{0,88},{4,88}},
          color={107,45,134}));
      connect(R0.port_a,oxyR1. products[1]) annotation (Line(
          points={{0,88},{0,74},{-6,74}},
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
          points={{-10,74},{-42,74},{-42,-36}},
          color={107,45,134}));
      connect(oxyR2.products[2],oxygen_unbound. port_a)
                                          annotation (Line(
          points={{-12,32},{-42,32},{-42,-36}},
          color={107,45,134}));
      connect(oxyR3.products[2],oxygen_unbound. port_a)
                                          annotation (Line(
          points={{-12,-14},{-42,-14},{-42,-36}},
          color={107,45,134}));
      connect(oxyR4.products[2],oxygen_unbound. port_a)
                                          annotation (Line(
          points={{-12,-56},{-42,-56},{-42,-36}},
          color={107,45,134}));
      connect(oxygen_unbound.port_a, oxyT1.products[2])
                                          annotation (Line(
          points={{-42,-36},{-42,74},{46,74}},
          color={107,45,134}));
      connect(oxygen_unbound.port_a, oxyT2.products[2])
                                          annotation (Line(
          points={{-42,-36},{-42,32},{46,32}},
          color={107,45,134}));
      connect(oxygen_unbound.port_a, oxyT3.products[2])
                                          annotation (Line(
          points={{-42,-36},{-42,-14},{46,-14}},
          color={107,45,134}));
      connect(oxygen_unbound.port_a, oxyT4.products[2])
                                          annotation (Line(
          points={{-42,-36},{-42,-56},{46,-56}},
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
          points={{54,0},{54,-14},{42,-14}},
          color={158,66,200}));
      connect(T2.solution, solution.solution) annotation (Line(points={{38,-10},
              {38,-99.74},{60.8,-99.74}},
                                color={127,127,0}));
      connect(substance.solution, solution.solution) annotation (Line(points={{
              70,-72},{66,-72},{66,-99.74},{60.8,-99.74}}, color={127,127,0}));
      connect(clock.y, O2_in_air.partialPressure)
        annotation (Line(points={{-84,51},{-84,32},{-84,32}}, color={0,0,127}));
      annotation ( Documentation(info="<html>
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
<p><i>2013-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Allosteric_Hemoglobin_MWC;

    model Allosteric_Hemoglobin2_MWC
      "Monod,Wyman,Changeux (1965) - The same allosteric hemoglobin model as Allosteric_Hemoglobin_MWC implemented by Speciation blocks"
      extends Modelica.Icons.Example;

      constant Modelica.Units.SI.AmountOfSubstance THb=0.001
        "Total amount of hemoglobin";

      constant Modelica.Units.SI.Temperature T=298.15 "Base Temperature";
      constant Real RT=Modelica.Constants.R*T;

      // constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
      //   "Amount of solution used for molarity to mole fraction conversion";
      constant Modelica.Units.SI.Volume OneLiter=0.001;

      parameter Real L=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.Units.SI.Concentration KR=0.000671946*(55.508/
          38.7) "oxygen dissociation on relaxed(R) hemoglobin subunit";
                                                                  //*7.875647668393782383419689119171e-5
      //10.500001495896 7.8756465463794e-05
      parameter Modelica.Units.SI.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.Units.SI.MoleFraction KRx=KR*OneLiter;
      //AmountOfSolutionIn1L;
      parameter Modelica.Units.SI.MoleFraction KTx=KT*OneLiter;
      //AmountOfSolutionIn1L;
      parameter Modelica.Units.SI.ChemicalPotential DfG_O2=-RT*log(0.0013);
      parameter Modelica.Units.SI.ChemicalPotential DfG_uR=0;
      parameter Modelica.Units.SI.ChemicalPotential DfG_uRO2=DfG_uR +
          DfG_O2 + RT*log(KRx);
      parameter Modelica.Units.SI.ChemicalPotential DfG_uT=0;
      parameter Modelica.Units.SI.ChemicalPotential DfG_uTO2=DfG_uT +
          DfG_O2 + RT*log(KTx);
      parameter Modelica.Units.SI.ChemicalPotential DfG_tT=0;
      parameter Modelica.Units.SI.ChemicalPotential DfG_tR=DfG_tT + RT*
          log(L);

      parameter Real KC = 1e-6 "Slow down factor";
                               //0.000001
      Chemical.Components.Solution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,42}})));

      Chemical.Components.Reaction quaternaryForm(
        nS=1,
        nP=1,
        KC=KC)
        annotation (Placement(transformation(extent={{12,-58},{32,-38}})));
      Chemical.Components.Speciation R0_in_R(NumberOfSubunits=4)
        annotation (Placement(transformation(extent={{-46,-48},{-26,-28}})));
       // AmountOfSubstance_start=4e-11)
      Chemical.Components.Speciation T0_in_T(NumberOfSubunits=4)
        annotation (Placement(transformation(extent={{76,-48},{56,-28}})));
       // AmountOfSubstance_start=totalAmountOfHemoglobin)
      Chemical.Components.Substance OxyRHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_O2 + RT*log(KRx) + DfG_tR/4),
        each use_mass_start=false,
        each amountOfSubstance_start=5.88e-9)
        "Oxygenated subunit in R structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{-96,-18},{-76,2}})));

      Chemical.Components.Reaction oxygenation_R[4](
        each nS=1,
        each nP=2,
        each KC=KC)
        annotation (Placement(transformation(extent={{-68,-18},{-48,2}})));
      Chemical.Components.Substance DeoxyRHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_tR/4),
        each use_mass_start=false,
        each amountOfSubstance_start=1.58e-7)
        "Deoxygenated subunit in R structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{-8,-20},{-28,0}})));

      Chemical.Components.Substance OxyTHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_O2 + RT*log(KTx) + DfG_tT/4),
        each use_mass_start=false,
        each amountOfSubstance_start=1e-4)
        "Oxygenated subunit in T structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{14,-18},{34,2}})));

      Chemical.Components.Reaction oxygenation_T[4](
        each nS=1,
        each nP=2,
        each KC=KC)
        annotation (Placement(transformation(extent={{42,-18},{62,2}})));
      Chemical.Components.Substance DeoxyTHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_tT/4),
        each use_mass_start=false,
        each amountOfSubstance_start=THb - 1e-4 - 1.58e-7 - 5.88e-9)
        "Deoxygenated subunit in T structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{96,-20},{76,0}})));

      Chemical.Components.Substance oxygen_unbound(
        substanceData(DfG=DfG_O2),
        use_mass_start=false,
        amountOfSubstance_start=2e-9)
        annotation (Placement(transformation(extent={{-2,6},{18,26}})));
      Modelica.Blocks.Sources.ContinuousClock clock(offset=1) annotation (
         Placement(transformation(extent={{-40,74},{-20,94}})));
      Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
          usePartialPressureInput=true, substanceData=
            Chemical.Substances.Oxygen_gas()) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={8,68})));
      Chemical.Components.GasSolubility partialPressure1(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={8,40})));

      Real sO2 "Hemoglobin oxygen saturation";
      Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
          mass_start=1)
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
          points={{62,-6},{70,-10},{76,-10}},
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
          points={{-48,-6},{-38,-6},{-38,-10},{-28,-10}},
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
          points={{62,-10},{70,-10},{70,16},{18,16}},
          color={107,45,134}));
        connect(oxygenation_R[i].products[2], oxygen_unbound.port_a) annotation (Line(
          points={{-48,-10},{-34,-10},{-34,16},{18,16}},
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
<p><i>2013-2018</i></p>
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

      parameter Modelica.Units.SI.MolarEnthalpy Ho=59000
        "Enthalpy of deoxygenation";
      parameter Modelica.Units.SI.MoleFraction Ko37
        "KRx and KTx at 37degC";
      parameter Modelica.Units.SI.MoleFraction Ko25=Ko37*exp((Ho/Modelica.Constants.R)
          *(1/310.15 - 1/298.15)) "KRx and KTx at 25degC";

      parameter Modelica.Units.SI.MolarEnthalpy Hco=59000
        "Enthalpy of carbon monoxide dissociation";
      parameter Modelica.Units.SI.MoleFraction Kco37
        "Carboxyhemoglobin dissociation at 37degC";
      parameter Modelica.Units.SI.MoleFraction Kco25=Kco37*exp((Hco/
          Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "Carboxyhemoglobin dissociation at 25degC";

      parameter Modelica.Units.SI.MolarEnthalpy Hh
        "Enthalpy of deprotonation of h site";
      parameter Modelica.Units.SI.MoleFraction Kh37
        "KRhx and KThx at 37 degC";
      parameter Modelica.Units.SI.MoleFraction Kh25=Kh37*exp(((Hh)/
          Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRhx and KThx at 25 degC";

      parameter Modelica.Units.SI.MolarEnthalpy Hz
        "Enthalpy of deprotonation of -NH3+ terminus";
      parameter Modelica.Units.SI.MoleFraction Kz37
        "KRzx and KTzx at 37 degC";
      parameter Modelica.Units.SI.MoleFraction Kz25=Kz37*exp(((Hz)/
          Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRzx and KTzx at 25 degC";

      parameter Modelica.Units.SI.MolarEnthalpy Hc
        "Enthalpy of carboxylation";
      parameter Modelica.Units.SI.MoleFraction Kc37
        "KRcx and KTcx at 37degC";
      parameter Modelica.Units.SI.MoleFraction Kc25=Kc37*exp((Hc/Modelica.Constants.R)
          *(1/310.15 - 1/298.15)) "KRcx and KTcx at 25degC";

      parameter Modelica.Units.SI.ChemicalPotential DfG_O2=-RT*log(0.0013)
           + 0;
      parameter Modelica.Units.SI.ChemicalPotential DfH_O2=0;

      parameter Modelica.Units.SI.ChemicalPotential DfG_CO=-RT*log(
          0.00099) - 137300;                                                                       //==Chemical.Examples.Substances.CarbonMonoxide_aqueous.DfG
      parameter Modelica.Units.SI.ChemicalPotential DfH_CO=-276900;

      parameter Modelica.Units.SI.ChemicalPotential DfG_CO2=-RT*log(0.034)
           - 394400;
      parameter Modelica.Units.SI.ChemicalPotential DfH_CO2=-412900;

      parameter Modelica.Units.SI.ChemicalPotential DfG_selectedForm
        "DfG_tR and DfG_tT";
      parameter Modelica.Units.SI.MolarEnthalpy DfH_selectedForm=0
        "DfH_tR and DfH_tT";

      parameter Real KC = 1e-3 "Slow down factor";
                               //0.000001
      parameter Modelica.Units.SI.MoleFraction initialO2
        "Initial mole fraction of unbound oxygen disoluted around hemoglobin";
      parameter Modelica.Units.SI.MoleFraction initialH
        "Initial mole fraction of H+";
      parameter Modelica.Units.SI.MoleFraction initialCO2
        "Initial mole fraction of unbound carbon dioxide disoluted around hemoglobin";
      parameter Modelica.Units.SI.AmountOfSubstance initialHb
        "Initial amount of hemoglobin tetramers in this quaternary form";

      Chemical.Components.Speciation speciation(NumberOfSubunits=N)
        annotation (Placement(transformation(extent={{-18,-72},{2,-52}})));
       // AmountOfSubstance_start=4e-11)
      // AmountOfSubstance_start=totalAmountOfHemoglobin)
      Chemical.Components.Substance OxyHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_O2 + RT*log(Ko25) + DfG_selectedForm/N, DfH=DfH_O2
               - Ho + DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=initialO2*initialHb/(Ko37 + initialO2))
        "Oxygenated subunit"
        annotation (Placement(transformation(extent={{-88,14},{-68,34}})));

      Chemical.Components.Reaction o[4](
        each nS=1,
        each nP=2,
        each KC=KC)
        annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
      Chemical.Components.Substance DeoxyHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=Ko37*initialHb/(Ko37 + initialO2))
        "Deoxygenated subunit"
        annotation (Placement(transformation(extent={{-8,12},{-28,32}})));

      Chemical.Interfaces.SolutionPort solution annotation (Placement(
            transformation(extent={{-50,-82},{-32,-62}}), iconTransformation(
              extent={{-50,-90},{-30,-70}})));
      Chemical.Interfaces.SubstancePort_b O2 annotation (Placement(
            transformation(extent={{-28,32},{-8,52}}), iconTransformation(
              extent={{-90,70},{-70,90}})));
      Chemical.Interfaces.SubstancePort_a selectedForm annotation (Placement(
            transformation(extent={{26,-82},{46,-62}}), iconTransformation(
              extent={{30,-90},{50,-70}})));
      Chemical.Components.Substance HmAH[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=RT*log(Kh25) + DfG_selectedForm/N, DfH=-Hh +
              DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=initialH*initialHb/(Kh37 + initialH))
        "Protonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{10,12},{30,32}})));

      Chemical.Components.Reaction h[4](
        each nS=1,
        each nP=2,
        each KC=KC)
        annotation (Placement(transformation(extent={{36,32},{56,12}})));
      Chemical.Components.Substance HmA[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=Kh37*initialHb/(Kh37 + initialH))
        "Deprotonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{92,14},{72,34}})));

      Chemical.Components.Substance HmNH3[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=RT*log(Kz25) + DfG_selectedForm/N, DfH=-Hz +
              DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=(initialH^2)*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Protonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{-84,-42},{-64,-22}})));

      Chemical.Components.Reaction z[4](
        each nS=1,
        each nP=2,
        each KC=KC)
        annotation (Placement(transformation(extent={{-54,-42},{-34,-22}})));
      Chemical.Components.Substance HmNH2[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=initialH*Kz37*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Deprotonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{12,-44},{-8,-24}})));

      Chemical.Components.Reaction c[4](
        each nP=2,
        each KC=KC,
        each nS=2)
        annotation (Placement(transformation(extent={{20,-42},{40,-22}})));
      Chemical.Interfaces.SubstancePort_b CO2 annotation (Placement(
            transformation(extent={{-8,-26},{12,-6}}), iconTransformation(
              extent={{10,70},{30,90}})));
      Chemical.Components.Substance HmNHCOO[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_CO2 - RT*log(Kc25) + DfG_selectedForm/N, DfH=
              DfH_CO2 + Hc + DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=(Kz37*Kc37*initialCO2)*initialHb/(initialH^2 +
            initialH*Kz37 + Kz37*Kc37*initialCO2))
        "Carboxylated c site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{70,-44},{50,-24}})));

      Interfaces.SubstancePort_b H annotation (Placement(transformation(extent={{50,-12},
                {70,8}}),      iconTransformation(extent={{70,70},{90,90}})));
    equation

      connect(OxyHm.port_a, o.substrates[1])
        annotation (Line(points={{-68,24},{-60,24},{-60,24}}, color={158,66,200}));
      connect(o.products[1], DeoxyHm.port_a)
        annotation (Line(points={{-40,26},{-28,26},{-28,22}},
                                                     color={158,66,200}));

      for i in 1:4 loop
        connect(h[i].products[2], H) annotation (Line(
            points={{56,24},{60,24},{60,-2}},
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

        connect(o[i].products[2], O2) annotation (Line(points={{-40,22},{-16,22},{-16,
                42},{-18,42}}, color={158,66,200}));
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
            points={{-34,-34},{-22,-34},{-22,-2},{60,-2}},
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
            points={{40,-34},{46,-34},{46,-2},{60,-2}},
            color={158,66,200}));
        connect(CO2, c[i].substrates[2]) annotation (Line(
            points={{2,-16},{16,-16},{16,-34},{20,-34}},
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
          points={{56,20},{64,24},{72,24}},
          color={158,66,200}));

      connect(z.products[1], HmNH2.port_a) annotation (Line(
          points={{-34,-30},{-22,-30},{-22,-34},{-8,-34}},
          color={107,45,134}));

      connect(HmNH3.port_a, z.substrates[1]) annotation (Line(
          points={{-64,-32},{-62,-32},{-60,-32},{-54,-32}},
          color={158,66,200}));

      connect(HmNH2.port_a, c.substrates[1]) annotation (Line(
          points={{-8,-34},{6,-34},{6,-30},{20,-30}},
          color={158,66,200}));

      connect(HmNHCOO.port_a, c.products[1]) annotation (Line(
          points={{50,-34},{46,-34},{46,-30},{40,-30}},
          color={158,66,200}));

      connect(solution, solution) annotation (Line(
          points={{-41,-72},{-41,-72}},
          color={127,127,0}));
      connect(H, H) annotation (Line(
          points={{60,-2},{60,-2}},
          color={158,66,200}));

      annotation (
        Documentation(revisions="<html>
<p><i>2013-2018</i></p>
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

      constant Modelica.Units.SI.AmountOfSubstance THb=0.001
        "Total amount of hemoglobin";

      constant Real RT=Modelica.Constants.R*298.15;

      // constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
      //   "Amount of solution used for molarity to mole fraction conversion";
      constant Modelica.Units.SI.Volume OneLiter=0.001;

      parameter Modelica.Units.SI.Pressure pCO2=5330 "partial pressure of CO2";
      parameter Real pH=7.2 "initial pH";

      parameter Real L_old=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.Units.SI.Concentration KR=0.000671946*(55.508/
          38.7) "oxygen dissociation on relaxed(R) hemoglobin subunit";

      parameter Modelica.Units.SI.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.Units.SI.MoleFraction KRo37=KR*OneLiter;
      parameter Modelica.Units.SI.MoleFraction KTo37=KT*OneLiter;

      parameter Modelica.Units.SI.ChemicalPotential DfG_O2=-RT*log(0.0013);
      parameter Modelica.Units.SI.ChemicalPotential DfG_CO2=-RT*log(0.034)
           - 394400;

      parameter Modelica.Units.SI.ChemicalPotential DfG_tT=0;
      parameter Modelica.Units.SI.ChemicalPotential DfG_tR=DfG_tT + RT*
          log(L);

      parameter Real KC = 1e-8 "Slow down factor";

      parameter Modelica.Units.SI.MoleFraction initialO2=1.9594e-07
        "Initial O2 at 37degC, pO2=100Pa";                         //at 25degC: 2.342e-8;
      parameter Modelica.Units.SI.MoleFraction initialH=10^(-pH);
      parameter Modelica.Units.SI.MoleFraction initialCO2=0.00136212*(pCO2/5330)
        "Initial CO2 at 37degC, pCO2=40mmHg";                                   //2.4217e-10
                                                                   //at 25degC: 3.267e-5;
      parameter Modelica.Units.SI.MoleFraction initialCO=1e-11
        "Initial CO at 37degC, pCO=0mmHg";
      //at 25degC: 3.267e-5;
      parameter Modelica.Units.SI.MoleFraction KRh37=10^(-6.89);
      parameter Modelica.Units.SI.MoleFraction KTh37=10^(-7.52);

      parameter Modelica.Units.SI.MoleFraction KRz37=10^(-7.25);
      parameter Modelica.Units.SI.MoleFraction KTz37=10^(-7.73);

      parameter Modelica.Units.SI.MoleFraction KRc37=(10^(-8.35))/(
          OneLiter);
      parameter Modelica.Units.SI.MoleFraction KTc37=(10^(-7.54))/(
          OneLiter);

      parameter Real L=L_old
        *
        (((KTh37/((10^(-7.2))+KTh37)) / (KRh37/((10^(-7.2))+KRh37)))^4)
        *
        (((KTz37*((10^(-7.2))^2 + KRz37*(10^(-7.2)) + KRz37*KRc37*(2.4217e-5)))/(KRz37*((10^(-7.2))^2 + KTz37*(10^(-7.2)) + KTz37*KTc37*(2.4217e-5))))^4)
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";

      Chemical.Components.Solution solution(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-56},{100,32}})));

      Chemical.Components.Reaction quaternaryForm(
        nS=1,
        nP=1,
        KC=KC)
        annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

      Chemical.Components.Substance O2_free(
        substanceData(DfG=DfG_O2, DfH=-11700),
        use_mass_start=false,
        amountOfSubstance_start=initialO2)
        annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
      Modelica.Blocks.Sources.ContinuousClock oxygenSource(offset=1000)
        annotation (Placement(transformation(extent={{-78,48},{-58,68}})));
      Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
        usePartialPressureInput=true,
        substanceData=Chemical.Substances.Oxygen_gas(),
        Temperature=310.15) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={-36,58})));
      Chemical.Components.GasSolubility partialPressure1(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,32})));

      Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
          mass_start=1)
        annotation (Placement(transformation(extent={{56,-54},{76,-34}})));
    //   (AmountOfSolutionIn1L -  THb - (initialO2 + initialCO2)*AmountOfSolutionIn1L)/55.508)
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
      Chemical.Sources.ExternalMoleFraction H(
        substanceData=Chemical.Substances.Proton_aqueous(),
        MoleFraction=initialH,
        Temperature=310.15) annotation (Placement(transformation(extent={{10,-10},
                {-10,10}}, origin={-12,-18})));
      Chemical.Components.Substance CO2_free(
        substanceData(DfG=DfG_CO2, DfH=-412900),
        use_mass_start=false,
        amountOfSubstance_start=initialCO2)
        annotation (Placement(transformation(extent={{86,-8},{66,12}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
        substanceData=Chemical.Substances.CarbonDioxide_gas(),
        PartialPressure(displayUnit="Pa") = pCO2,
        Temperature=310.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={62,60})));
      Chemical.Components.GasSolubility partialPressure2(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,32})));

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
<p><i>2013-2018</i></p>
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
                100}})));
    end HemoglobinMultipleAllostery;

    model HemoglobinQuaternaryFormCO
      "Hemoglobib quaternary form - part of multiple-ligand allosteric hemoglobin model"

      constant Integer N=12
        "Number of distinguished independent sides in quaternary structure";
      constant Real RT=Modelica.Constants.R*298.15;

      parameter Modelica.Units.SI.MolarEnthalpy Ho=59000
        "Enthalpy of deoxygenation";
      parameter Modelica.Units.SI.MoleFraction Ko37
        "KRx and KTx at 37degC";
      parameter Modelica.Units.SI.MoleFraction Ko25=Ko37*exp((Ho/Modelica.Constants.R)
          *(1/310.15 - 1/298.15)) "KRx and KTx at 25degC";

      parameter Modelica.Units.SI.MolarEnthalpy Hco=59000
        "Enthalpy of carbon monoxide dissociation";
      parameter Modelica.Units.SI.MoleFraction Kco37
        "Carboxyhemoglobin dissociation at 37degC";
      parameter Modelica.Units.SI.MoleFraction Kco25=Kco37*exp((Hco/
          Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "Carboxyhemoglobin dissociation at 25degC";

      parameter Modelica.Units.SI.MolarEnthalpy Hh
        "Enthalpy of deprotonation of h site";
      parameter Modelica.Units.SI.MoleFraction Kh37
        "KRhx and KThx at 37 degC";
      parameter Modelica.Units.SI.MoleFraction Kh25=Kh37*exp(((Hh)/
          Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRhx and KThx at 25 degC";

      parameter Modelica.Units.SI.MolarEnthalpy Hz
        "Enthalpy of deprotonation of -NH3+ terminus";
      parameter Modelica.Units.SI.MoleFraction Kz37
        "KRzx and KTzx at 37 degC";
      parameter Modelica.Units.SI.MoleFraction Kz25=Kz37*exp(((Hz)/
          Modelica.Constants.R)*(1/310.15 - 1/298.15))
        "KRzx and KTzx at 25 degC";

      parameter Modelica.Units.SI.MolarEnthalpy Hc
        "Enthalpy of carboxylation";
      parameter Modelica.Units.SI.MoleFraction Kc37
        "KRcx and KTcx at 37degC";
      parameter Modelica.Units.SI.MoleFraction Kc25=Kc37*exp((Hc/Modelica.Constants.R)
          *(1/310.15 - 1/298.15)) "KRcx and KTcx at 25degC";

      parameter Modelica.Units.SI.ChemicalPotential DfG_O2=-RT*log(0.0013)
           + 0;
      parameter Modelica.Units.SI.ChemicalPotential DfH_O2=0;

      parameter Modelica.Units.SI.ChemicalPotential DfG_CO=-RT*log(
          0.00099) - 137300;                                                                       //==Chemical.Examples.Substances.CarbonMonoxide_aqueous.DfG
      parameter Modelica.Units.SI.ChemicalPotential DfH_CO=-276900;

      parameter Modelica.Units.SI.ChemicalPotential DfG_CO2=-RT*log(0.034)
           - 394400;
      parameter Modelica.Units.SI.ChemicalPotential DfH_CO2=-412900;

      parameter Modelica.Units.SI.ChemicalPotential DfG_selectedForm
        "DfG_tR and DfG_tT";
      parameter Modelica.Units.SI.MolarEnthalpy DfH_selectedForm=0
        "DfH_tR and DfH_tT";

      parameter Real KC = 1e-3 "Slow down factor";
                               //0.000001
      parameter Modelica.Units.SI.MoleFraction initialO2
        "Initial mole fraction of unbound oxygen disoluted around hemoglobin";
      parameter Modelica.Units.SI.MoleFraction initialCO
        "Initial mole fraction of unbound carbon monoxide disoluted around hemoglobin";
      parameter Modelica.Units.SI.MoleFraction initialH
        "Initial mole fraction of H+";
      parameter Modelica.Units.SI.MoleFraction initialCO2
        "Initial mole fraction of unbound carbon dioxide disoluted around hemoglobin";
      parameter Modelica.Units.SI.AmountOfSubstance initialHb
        "Initial amount of hemoglobin tetramers in this quaternary form";

      Chemical.Components.Speciation speciation(NumberOfSubunits=N)
        annotation (Placement(transformation(extent={{-18,-72},{2,-52}})));
       // AmountOfSubstance_start=4e-11)
      // AmountOfSubstance_start=totalAmountOfHemoglobin)
      Chemical.Components.Substance OxyHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_O2 + RT*log(Ko25) + DfG_selectedForm/N, DfH=DfH_O2
               - Ho + DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=(initialO2/Ko37)*initialHb/(1 + initialO2/Ko37 +
            initialCO/Kco37)) "Oxygenated subunit"
        annotation (Placement(transformation(extent={{-88,14},{-68,34}})));

      Chemical.Components.Reaction o[4](
        each nS=1,
        each nP=2,
        each KC=KC)
        annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
      Chemical.Components.Substance DeoxyHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=initialHb/(1 + initialO2/Ko37 + initialCO/Kco37))
        "Deoxygenated subunit"
        annotation (Placement(transformation(extent={{-8,12},{-28,32}})));

      Chemical.Interfaces.SolutionPort solution annotation (Placement(
            transformation(extent={{-50,-82},{-32,-62}}), iconTransformation(
              extent={{-50,-90},{-30,-70}})));
      Chemical.Interfaces.SubstancePort_b O2 annotation (Placement(
            transformation(extent={{-28,32},{-8,52}}), iconTransformation(
              extent={{-90,70},{-70,90}})));
      Chemical.Interfaces.SubstancePort_a selectedForm annotation (Placement(
            transformation(extent={{26,-82},{46,-62}}), iconTransformation(
              extent={{30,-90},{50,-70}})));
      Chemical.Components.Substance HmAH[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=RT*log(Kh25) + DfG_selectedForm/N, DfH=-Hh +
              DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=initialH*initialHb/(Kh37 + initialH))
        "Protonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{10,12},{30,32}})));

      Chemical.Components.Reaction h[4](
        each nS=1,
        each nP=2,
        each KC=KC)
        annotation (Placement(transformation(extent={{36,32},{56,12}})));
      Chemical.Components.Substance HmA[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=Kh37*initialHb/(Kh37 + initialH))
        "Deprotonated h site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{92,14},{72,34}})));

      Chemical.Components.Substance HmNH3[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=RT*log(Kz25) + DfG_selectedForm/N, DfH=-Hz +
              DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=(initialH^2)*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Protonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{-84,-42},{-64,-22}})));

      Chemical.Components.Reaction z[4](
        each nS=1,
        each nP=2,
        each KC=KC)
        annotation (Placement(transformation(extent={{-54,-42},{-34,-22}})));
      Chemical.Components.Substance HmNH2[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=initialH*Kz37*initialHb/(initialH^2 + initialH*
            Kz37 + Kz37*Kc37*initialCO2))
        "Deprotonated z site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{12,-44},{-8,-24}})));

      Chemical.Components.Reaction c[4](
        each nP=2,
        each KC=KC,
        each nS=2)
        annotation (Placement(transformation(extent={{20,-42},{40,-22}})));
      Chemical.Interfaces.SubstancePort_b CO2 annotation (Placement(
            transformation(extent={{-8,-26},{12,-6}}), iconTransformation(
              extent={{10,70},{30,90}})));
      Chemical.Components.Substance HmNHCOO[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfG=DfG_CO2 - RT*log(Kc25) + DfG_selectedForm/N, DfH=
              DfH_CO2 + Hc + DfH_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=(Kz37*Kc37*initialCO2)*initialHb/(initialH^2 +
            initialH*Kz37 + Kz37*Kc37*initialCO2))
        "Carboxylated c site of subunit in quaternary structure of hemoglobin tetramer"
        annotation (Placement(transformation(extent={{70,-44},{50,-24}})));

      Interfaces.SubstancePort_b H annotation (Placement(transformation(extent={{50,-12},
                {70,8}}),      iconTransformation(extent={{70,70},{90,90}})));
      Chemical.Components.Substance COHm[4](
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        each substanceData(DfH=DfH_CO - Hco + DfH_selectedForm/N, DfG=DfG_CO + RT*log(
              Kco25) + DfG_selectedForm/N),
        each use_mass_start=false,
        each amountOfSubstance_start=(initialCO/Kco37)*initialHb/(1 + initialO2/Ko37
             + initialCO/Kco37)) "Subunit with Carbon Monoxide"
        annotation (Placement(transformation(extent={{78,44},{58,64}})));

      Chemical.Components.Reaction o1[4](
        each nS=1,
        each nP=2,
        each KC=KC)
        annotation (Placement(transformation(extent={{46,44},{26,64}})));
      Chemical.Interfaces.SubstancePort_b CO annotation (Placement(
            transformation(extent={{-12,62},{8,82}}), iconTransformation(
              extent={{-50,70},{-30,90}})));
    equation

      connect(OxyHm.port_a, o.substrates[1])
        annotation (Line(points={{-68,24},{-60,24}},          color={158,66,200}));
      connect(o.products[1], DeoxyHm.port_a)
        annotation (Line(points={{-40,26},{-34,26},{-34,22},{-28,22}},
                                                     color={158,66,200}));

      for i in 1:4 loop
        connect(h[i].products[2], H) annotation (Line(
            points={{56,24},{60,24},{60,-2}},
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

        connect(o[i].products[2], O2) annotation (Line(points={{-40,22},{-16,22},{-16,
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
            points={{-34,-34},{-22,-34},{-22,-2},{60,-2}},
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
            points={{40,-34},{46,-34},{46,-2},{60,-2}},
            color={158,66,200}));
        connect(CO2, c[i].substrates[2]) annotation (Line(
            points={{2,-16},{16,-16},{16,-34},{20,-34}},
            color={158,66,200}));
        connect(HmNHCOO[i].solution, speciation.subunitSolution) annotation (Line(
          points={{66,-44},{-4,-44},{-4,-56}},
          color={127,127,0}));

        connect(COHm[i].solution, speciation.subunitSolution) annotation (Line(
          points={{74,44},{74,42},{94,42},{94,-44},{-4,-44},{-4,-56}},
          color={127,127,0},
          smooth=Smooth.None));
        connect(o1[i].products[2], CO) annotation (Line(
          points={{26,52},{14,52},{14,72},{-2,72}},
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
          points={{56,20},{64,24},{72,24}},
          color={158,66,200}));

      connect(z.products[1], HmNH2.port_a) annotation (Line(
          points={{-34,-30},{-22,-30},{-22,-34},{-8,-34}},
          color={107,45,134}));

      connect(HmNH3.port_a, z.substrates[1]) annotation (Line(
          points={{-64,-32},{-62,-32},{-60,-32},{-54,-32}},
          color={158,66,200}));

      connect(HmNH2.port_a, c.substrates[1]) annotation (Line(
          points={{-8,-34},{6,-34},{6,-30},{20,-30}},
          color={158,66,200}));

      connect(HmNHCOO.port_a, c.products[1]) annotation (Line(
          points={{50,-34},{46,-34},{46,-30},{40,-30}},
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
          points={{-28,22},{-2,22},{-2,56},{26,56}},
          color={158,66,200},
          smooth=Smooth.None));
      annotation (
        Documentation(revisions="<html>
<p><i>2013-2018</i></p>
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

      constant Modelica.Units.SI.AmountOfSubstance THb=0.001
        "Total amount of hemoglobin";

      constant Real RT=Modelica.Constants.R*298.15;

      // constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
      //   "Amount of solution used for molarity to mole fraction conversion";
      constant Modelica.Units.SI.Volume OneLiter=0.001;

      parameter Real L_old=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.Units.SI.Concentration KR=0.000671946*(55.508/
          38.7) "oxygen dissociation on relaxed(R) hemoglobin subunit";

      parameter Modelica.Units.SI.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.Units.SI.MoleFraction KRo37=KR*OneLiter;
      parameter Modelica.Units.SI.MoleFraction KTo37=KT*OneLiter;

      parameter Modelica.Units.SI.ChemicalPotential DfG_O2=-RT*log(0.0013);
      parameter Modelica.Units.SI.ChemicalPotential DfG_CO2=-RT*log(0.034)
           - 394400;

      parameter Modelica.Units.SI.ChemicalPotential DfG_tT=0;
      parameter Modelica.Units.SI.ChemicalPotential DfG_tR=DfG_tT + RT*
          log(L);

      parameter Real KC = 1e-8 "Slow down factor";

      parameter Modelica.Units.SI.MoleFraction initialO2=1.9594e-07
        "Initial O2 at 37degC, pO2=100Pa";                         //at 25degC: 2.342e-8;
      parameter Modelica.Units.SI.MoleFraction initialH=10^(-7.2);
      parameter Modelica.Units.SI.MoleFraction initialCO2=2.4217e-10
        "Initial CO2 at 37degC, pCO2=40mmHg";                      //at 25degC: 3.267e-5;
      parameter Modelica.Units.SI.MoleFraction initialCO=1e-12
        "Initial CO at 37degC, pCO=0mmHg";
      //at 25degC: 3.267e-5;
      parameter Modelica.Units.SI.MoleFraction KRh37=10^(-6.89);
      parameter Modelica.Units.SI.MoleFraction KTh37=10^(-7.52);

      parameter Modelica.Units.SI.MoleFraction KRz37=10^(-7.25);
      parameter Modelica.Units.SI.MoleFraction KTz37=10^(-7.73);

      parameter Modelica.Units.SI.MoleFraction KRc37=(10^(-8.35))/(
          OneLiter);
      parameter Modelica.Units.SI.MoleFraction KTc37=(10^(-7.54))/(
          OneLiter);

      parameter Real L=L_old
        *
        (((KTh37/((10^(-7.2))+KTh37)) / (KRh37/((10^(-7.2))+KRh37)))^4)
        *
        (((KTz37*((10^(-7.2))^2 + KRz37*(10^(-7.2)) + KRz37*KRc37*(2.4217e-5)))/(KRz37*((10^(-7.2))^2 + KTz37*(10^(-7.2)) + KTz37*KTc37*(2.4217e-5))))^4)
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";

      Chemical.Components.Solution solution(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-56},{100,32}})));

      Chemical.Components.Reaction quaternaryForm(
        nS=1,
        nP=1,
        KC=KC)
        annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

      Chemical.Components.Substance O2_free(
        substanceData(DfG=DfG_O2, DfH=-11700),
        use_mass_start=false,
        amountOfSubstance_start=initialO2)
        annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
      Modelica.Blocks.Sources.ContinuousClock oxygenSource(offset=2000)
        annotation (Placement(transformation(extent={{-78,48},{-58,68}})));
      Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
        usePartialPressureInput=true,
        substanceData=Chemical.Substances.Oxygen_gas(),
        Temperature=310.15) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={-36,58})));
      Chemical.Components.GasSolubility partialPressure1(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,32})));

      Components.Substance H2O(substanceData=Substances.Water_liquid(),
          mass_start=1)
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
      Chemical.Sources.ExternalMoleFraction H(
        substanceData=Chemical.Substances.Proton_aqueous(),
        MoleFraction=initialH,
        Temperature=310.15) annotation (Placement(transformation(extent={{10,-10},
                {-10,10}}, origin={-12,-18})));
      Chemical.Components.Substance CO2_free(
        substanceData(DfG=DfG_CO2, DfH=-412900),
        use_mass_start=false,
        amountOfSubstance_start=initialCO2)
        annotation (Placement(transformation(extent={{86,-8},{66,12}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
        substanceData=Chemical.Substances.CarbonDioxide_gas(),
        PartialPressure(displayUnit="kPa") = 5330,
        Temperature=310.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={62,60})));
      Chemical.Components.GasSolubility partialPressure2(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,32})));

      Real sCO "Hemoglobin carbon monoxide saturation";
      Real sO2 "Hemoglobin oxygen saturation";
      Real sCO2 "Hemoglobin carbon dioxide saturation";
      Real dH "Hemoglobin charge change caused by binding of Bohr's protons";
      Chemical.Components.Substance CO_free(
        substanceData=Chemical.Substances.CarbonMonoxide_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=initialCO)
        annotation (Placement(transformation(extent={{-92,8},{-72,28}})));
      Chemical.Components.GasSolubility partialPressure3(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={26,32})));
      Chemical.Sources.ExternalIdealGasSubstance CO_gas(
        substanceData=Chemical.Substances.CarbonMonoxide_gas(),
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
<p><i>2013-2018</i></p>
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

    model Joels57
      extends Modelica.Icons.Example;
       HemoglobinMultipleAllostery pCO2_2000(pH=7.6, pCO2(displayUnit="kPa") = 2000)
        annotation (Placement(transformation(extent={{-50,-18},{-30,2}})));
      HemoglobinMultipleAllostery pCO2_5330(pH=7.5, pCO2(displayUnit="kPa") = 5330)
        annotation (Placement(transformation(extent={{-6,-16},{14,4}})));
      HemoglobinMultipleAllostery pCO2_9330(pH=7.4, pCO2(displayUnit="kPa") = 9330)
        annotation (Placement(transformation(extent={{40,-16},{60,4}})));
      annotation (experiment(
          StopTime=70,
          __Dymola_NumberOfIntervals=50000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Dassl"));
    end Joels57;

    model HemoglobinTitration
      "Multiple-ligand allosteric hemoglobin model"
      extends Modelica.Icons.Example;

      constant Modelica.Units.SI.AmountOfSubstance THb=0.001
        "Total amount of hemoglobin";

      constant Real RT=Modelica.Constants.R*298.15;

      // constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
      //   "Amount of solution used for molarity to mole fraction conversion";
      constant Modelica.Units.SI.Volume OneLiter=0.001;

      parameter Real L_old=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.Units.SI.Concentration KR=0.000671946*(55.508/
          38.7) "oxygen dissociation on relaxed(R) hemoglobin subunit";

      parameter Modelica.Units.SI.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.Units.SI.MoleFraction KRo37=KR*OneLiter;
      parameter Modelica.Units.SI.MoleFraction KTo37=KT*OneLiter;

      parameter Modelica.Units.SI.ChemicalPotential DfG_O2=-RT*log(0.0013);
      parameter Modelica.Units.SI.ChemicalPotential DfG_CO2=-RT*log(0.034)
           - 394400;

      parameter Modelica.Units.SI.ChemicalPotential DfG_tT=0;
      parameter Modelica.Units.SI.ChemicalPotential DfG_tR=DfG_tT + RT*
          log(L);

      parameter Real KC = 1e-8 "Slow down factor";

      parameter Modelica.Units.SI.MoleFraction initialO2=0.000132233
        "Initial O2 at 37degC, pO2=100Pa";                         //at 25degC: 2.342e-8;
      parameter Modelica.Units.SI.MoleFraction initialH=10^(-6.9);
      parameter Modelica.Units.SI.MoleFraction initialCO2=0.00136212
        "Initial CO2 at 37degC, pCO2=40mmHg";
                                                                   //at 25degC: 3.267e-5;
      parameter Modelica.Units.SI.MoleFraction initialCO=1e-11
        "Initial CO at 37degC, pCO=0mmHg";
      //at 25degC: 3.267e-5;
      parameter Modelica.Units.SI.MoleFraction KRh37=10^(-6.89);
      parameter Modelica.Units.SI.MoleFraction KTh37=10^(-7.52);

      parameter Modelica.Units.SI.MoleFraction KRz37=10^(-7.25);
      parameter Modelica.Units.SI.MoleFraction KTz37=10^(-7.73);

      parameter Modelica.Units.SI.MoleFraction KRc37=(10^(-8.35))/(
          OneLiter);
      parameter Modelica.Units.SI.MoleFraction KTc37=(10^(-7.54))/(
          OneLiter);

      parameter Real L=L_old
        *
        (((KTh37/((10^(-7.2))+KTh37)) / (KRh37/((10^(-7.2))+KRh37)))^4)
        *
        (((KTz37*((10^(-7.2))^2 + KRz37*(10^(-7.2)) + KRz37*KRc37*(2.4217e-5)))/(KRz37*((10^(-7.2))^2 + KTz37*(10^(-7.2)) + KTz37*KTc37*(2.4217e-5))))^4)
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";

      Chemical.Components.Solution solution(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-96,-60},{104,28}})));

      Chemical.Components.Reaction quaternaryForm(
        nS=1,
        nP=1,
        KC=KC)
        annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

      Chemical.Components.Substance O2_free(
        substanceData(DfG=DfG_O2, DfH=-11700),
        use_mass_start=false,
        amountOfSubstance_start=initialO2)
        annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
      Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
        usePartialPressureInput=false,
        substanceData=Chemical.Substances.Oxygen_gas(),
        PartialPressure(displayUnit="mmHg") = 11999.01486735,
        Temperature=310.15) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={-36,58})));
      Chemical.Components.GasSolubility partialPressure1(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,32})));

      Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
          mass_start=1)
        annotation (Placement(transformation(extent={{56,-54},{76,-34}})));
    //   (AmountOfSolutionIn1L -  THb - (initialO2 + initialCO2)*AmountOfSolutionIn1L)/55.508)
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
      Chemical.Sources.ExternalMoleFraction H(
        substanceData=Chemical.Substances.Proton_aqueous(),
        MoleFraction=initialH,
        useMoleFractionInput=true,
        Temperature=310.15) annotation (Placement(transformation(extent={{10,-10},
                {-10,10}}, origin={-12,-18},
            rotation=270)));
      Chemical.Components.Substance CO2_free(
        substanceData(DfG=DfG_CO2, DfH=-412900),
        use_mass_start=false,
        amountOfSubstance_start=initialCO2)
        annotation (Placement(transformation(extent={{86,-8},{66,12}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
        substanceData=Chemical.Substances.CarbonDioxide_gas(),
        PartialPressure(displayUnit="mmHg") = 5332.8954966,
        Temperature=310.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={62,60})));
      Chemical.Components.GasSolubility partialPressure2(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,32})));

      Real sO2 "Hemoglobin oxygen saturation";
      Real sCO2 "Hemoglobin carbon dioxide saturation";
      Real dH "Hemoglobin charge change caused by binding of Bohr's protons";
      Modelica.Blocks.Sources.ContinuousClock pHSource(offset=690000)
        annotation (Placement(transformation(extent={{-54,-90},{-34,-70}})));
      Modelica.Blocks.Math.Exp exp1
        annotation (Placement(transformation(extent={{10,-90},{30,-70}})));
      Modelica.Blocks.Math.Gain gain(k=-log(10)/100000)
        annotation (Placement(transformation(extent={{-22,-90},{-2,-70}})));
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
            points={{-14,22},{-14,-2},{-56,-2}}, color={158,66,200}));

      connect(O2_free.solution, solution.solution) annotation (Line(points={{-72,-12},
              {-72,-54},{64,-54},{64,-59.12}},        color={127,127,0}));
      connect(solution.solution, H2O.solution) annotation (Line(
          points={{64,-59.12},{64,-54},{60,-54}},
          color={127,127,0}));

      connect(relaxed.solution, solution.solution) annotation (Line(
          points={{-48,-42},{-48,-54},{64,-54},{64,-59.12}},
          color={127,127,0}));
      connect(relaxed.O2, O2_free.port_a) annotation (Line(
          points={{-52,-26},{-52,-2},{-56,-2}},
          color={158,66,200}));
      connect(relaxed.selectedForm, quaternaryForm.substrates[1]) annotation (Line(
          points={{-40,-42},{-22,-42}},
          color={158,66,200}));
      connect(tensed.solution, solution.solution) annotation (Line(
          points={{26,-42},{26,-54},{64,-54},{64,-59.12}},
          color={127,127,0}));
      connect(tensed.O2, O2_free.port_a) annotation (Line(
          points={{30,-26},{30,-2},{-56,-2}},
          color={158,66,200}));
      connect(tensed.selectedForm, quaternaryForm.products[1]) annotation (Line(
          points={{18,-42},{-2,-42}},
          color={158,66,200}));
      connect(H.port_a, relaxed.H) annotation (Line(
          points={{-12,-8},{-36,-8},{-36,-26}},
          color={158,66,200}));
      connect(H.port_a, tensed.H) annotation (Line(
          points={{-12,-8},{14,-8},{14,-26}},
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
          points={{82,-8},{82,-59.12},{64,-59.12}},
          color={127,127,0}));
      connect(exp1.y, H.moleFractionInput) annotation (Line(points={{31,-80},{38,
              -80},{38,-34},{-12,-34},{-12,-28}}, color={0,0,127}));
      connect(pHSource.y,gain. u) annotation (Line(
          points={{-33,-80},{-24,-80}},
          color={0,0,127}));
      connect(exp1.u, gain.y)
        annotation (Line(points={{8,-80},{-1,-80}}, color={0,0,127}));
      annotation (          experiment(
          StopTime=100000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Dassl"),
        Documentation(revisions="<html>
<p><i>2013-2018</i></p>
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
                100}})));
    end HemoglobinTitration;

    model HemoglobinCarboxylation
      "Multiple-ligand allosteric hemoglobin model"
      extends Modelica.Icons.Example;

      constant Modelica.Units.SI.AmountOfSubstance THb=0.001
        "Total amount of hemoglobin";

      constant Real RT=Modelica.Constants.R*298.15;

      // constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
      //   "Amount of solution used for molarity to mole fraction conversion";
      constant Modelica.Units.SI.Volume OneLiter=0.001;

      parameter Real L_old=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
      parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
      parameter Modelica.Units.SI.Concentration KR=0.000671946*(55.508/
          38.7) "oxygen dissociation on relaxed(R) hemoglobin subunit";

      parameter Modelica.Units.SI.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

      parameter Modelica.Units.SI.MoleFraction KRo37=KR*OneLiter;
      parameter Modelica.Units.SI.MoleFraction KTo37=KT*OneLiter;

      parameter Modelica.Units.SI.ChemicalPotential DfG_O2=-RT*log(0.0013);
      parameter Modelica.Units.SI.ChemicalPotential DfG_CO2=-RT*log(0.034)
           - 394400;

      parameter Modelica.Units.SI.ChemicalPotential DfG_tT=0;
      parameter Modelica.Units.SI.ChemicalPotential DfG_tR=DfG_tT + RT*
          log(L);

      parameter Real KC = 1e-8 "Slow down factor";

      parameter Modelica.Units.SI.MoleFraction initialO2=0.000132233
        "Initial O2 at 37degC, pO2=100Pa";                         //at 25degC: 2.342e-8;
      parameter Modelica.Units.SI.MoleFraction initialH=10^(-6.9);
      parameter Modelica.Units.SI.MoleFraction initialCO2=0.00136212
        "Initial CO2 at 37degC, pCO2=40mmHg";
                                                                   //at 25degC: 3.267e-5;
      parameter Modelica.Units.SI.MoleFraction initialCO=1e-11
        "Initial CO at 37degC, pCO=0mmHg";
      //at 25degC: 3.267e-5;
      parameter Modelica.Units.SI.MoleFraction KRh37=10^(-6.89);
      parameter Modelica.Units.SI.MoleFraction KTh37=10^(-7.52);

      parameter Modelica.Units.SI.MoleFraction KRz37=10^(-7.25);
      parameter Modelica.Units.SI.MoleFraction KTz37=10^(-7.73);

      parameter Modelica.Units.SI.MoleFraction KRc37=(10^(-8.35))/(
          OneLiter);
      parameter Modelica.Units.SI.MoleFraction KTc37=(10^(-7.54))/(
          OneLiter);

      parameter Real L=L_old
        *
        (((KTh37/((10^(-7.2))+KTh37)) / (KRh37/((10^(-7.2))+KRh37)))^4)
        *
        (((KTz37*((10^(-7.2))^2 + KRz37*(10^(-7.2)) + KRz37*KRc37*(2.4217e-5)))/(KRz37*((10^(-7.2))^2 + KTz37*(10^(-7.2)) + KTz37*KTc37*(2.4217e-5))))^4)
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";

      Chemical.Components.Solution solution(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-96,-60},{104,28}})));

      Chemical.Components.Reaction quaternaryForm(
        nS=1,
        nP=1,
        KC=KC)
        annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

      Chemical.Components.Substance O2_free(
        substanceData(DfG=DfG_O2, DfH=-11700),
        use_mass_start=false,
        amountOfSubstance_start=initialO2)
        annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
      Chemical.Sources.ExternalIdealGasSubstance oxygen_in_air(
        usePartialPressureInput=false,
        substanceData=Chemical.Substances.Oxygen_gas(),
        PartialPressure(displayUnit="mmHg") = 11999.01486735,
        Temperature=310.15) annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, origin={-36,58})));
      Chemical.Components.GasSolubility partialPressure1(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,32})));

      Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
          mass_start=1)
        annotation (Placement(transformation(extent={{56,-54},{76,-34}})));
    //   (AmountOfSolutionIn1L -  THb - (initialO2 + initialCO2)*AmountOfSolutionIn1L)/55.508)
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
      Chemical.Sources.ExternalMoleFraction H(
        substanceData=Chemical.Substances.Proton_aqueous(),
        MoleFraction=initialH,
        useMoleFractionInput=false,
        Temperature=310.15) annotation (Placement(transformation(extent={{10,-10},
                {-10,10}}, origin={-12,-18},
            rotation=270)));
      Chemical.Components.Substance CO2_free(
        substanceData(DfG=DfG_CO2, DfH=-412900),
        use_mass_start=false,
        amountOfSubstance_start=initialCO2)
        annotation (Placement(transformation(extent={{86,-8},{66,12}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
        substanceData=Chemical.Substances.CarbonDioxide_gas(),
        usePartialPressureInput=true,
        PartialPressure(displayUnit="kPa") = 5330,
        Temperature=310.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={62,60})));
      Chemical.Components.GasSolubility partialPressure2(KC=KC) annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,32})));

      Real sO2 "Hemoglobin oxygen saturation";
      Real sCO2 "Hemoglobin carbon dioxide saturation";
      Real dH "Hemoglobin charge change caused by binding of Bohr's protons";
      Modelica.Blocks.Sources.ContinuousClock carbonDioxideSource(offset=
            1000)
        annotation (Placement(transformation(extent={{10,70},{30,90}})));
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
            points={{-14,22},{-14,-2},{-56,-2}}, color={158,66,200}));

      connect(O2_free.solution, solution.solution) annotation (Line(points={{-72,-12},
              {-72,-54},{64,-54},{64,-59.12}},        color={127,127,0}));
      connect(solution.solution, H2O.solution) annotation (Line(
          points={{64,-59.12},{64,-54},{60,-54}},
          color={127,127,0}));

      connect(relaxed.solution, solution.solution) annotation (Line(
          points={{-48,-42},{-48,-54},{64,-54},{64,-59.12}},
          color={127,127,0}));
      connect(relaxed.O2, O2_free.port_a) annotation (Line(
          points={{-52,-26},{-52,-2},{-56,-2}},
          color={158,66,200}));
      connect(relaxed.selectedForm, quaternaryForm.substrates[1]) annotation (Line(
          points={{-40,-42},{-22,-42}},
          color={158,66,200}));
      connect(tensed.solution, solution.solution) annotation (Line(
          points={{26,-42},{26,-54},{64,-54},{64,-59.12}},
          color={127,127,0}));
      connect(tensed.O2, O2_free.port_a) annotation (Line(
          points={{30,-26},{30,-2},{-56,-2}},
          color={158,66,200}));
      connect(tensed.selectedForm, quaternaryForm.products[1]) annotation (Line(
          points={{18,-42},{-2,-42}},
          color={158,66,200}));
      connect(H.port_a, relaxed.H) annotation (Line(
          points={{-12,-8},{-36,-8},{-36,-26}},
          color={158,66,200}));
      connect(H.port_a, tensed.H) annotation (Line(
          points={{-12,-8},{14,-8},{14,-26}},
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
          points={{82,-8},{82,-59.12},{64,-59.12}},
          color={127,127,0}));
      connect(carbonDioxideSource.y, CO2_gas.partialPressure) annotation (
          Line(points={{31,80},{62,80},{62,70}}, color={0,0,127}));
      annotation (          experiment(
          StopTime=100000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Dassl"),
        Documentation(revisions="<html>
<p><i>2013-2018</i></p>
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
                100}})));
    end HemoglobinCarboxylation;
  end Hemoglobin;

  package CheckSubstancesData
    model SimpleReaction
      "The simple chemical reaction A<->B with equilibrium B/A = 2"
       extends Modelica.Icons.Example;

      constant Real K = 2 "Dissociation constant of the reaction";

      constant Modelica.Units.SI.Temperature T_25degC=298.15
        "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

    Chemical.Sensors.DissociationCoefficient dissociationCoefficient(nS=1, nP=1)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    Chemical.Sources.PureSubstance A(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible, substanceData(DfG=0))
        annotation (Placement(transformation(extent={{-56,-10},{-36,10}})));
    Chemical.Sources.PureSubstance B(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible, substanceData(DfG=-R*T_25degC
              *log(K)))
        annotation (Placement(transformation(extent={{60,-10},{40,10}})));
      inner Modelica.Fluid.System system(p_ambient=100000, T_ambient=298.15)
        annotation (Placement(transformation(extent={{-58,62},{-38,82}})));
    equation
      connect(A.port_a, dissociationCoefficient.substrates[1])
        annotation (Line(points={{-36,0},{-10,0}}, color={158,66,200}));
      connect(dissociationCoefficient.products[1], B.port_a)
        annotation (Line(points={{10,0},{40,0}}, color={158,66,200}));
      annotation (Documentation(revisions="<html>
<p><i>2015-2018</i></p>
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

      constant Modelica.Units.SI.Temperature T_25degC=298.15
        "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

    Chemical.Sources.PureSubstance A
        annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
      Chemical.Sensors.DissociationCoefficient reaction(nS=2, nP=1)
        annotation (Placement(transformation(extent={{4,-8},{24,12}})));
    Chemical.Sources.PureSubstance B
        annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
    Chemical.Sources.PureSubstance C(substanceData(DfG=-R*T_25degC*log(Kx)))
        annotation (Placement(transformation(extent={{68,-8},{48,12}})));

    equation

      connect(A.port_a, reaction.substrates[1]) annotation (Line(points={
              {-14,12},{-4,12},{-4,1.5},{4,1.5}}, color={158,66,200}));
      connect(B.port_a, reaction.substrates[2]) annotation (Line(points={
              {-14,-14},{-4,-14},{-4,2.5},{4,2.5}}, color={158,66,200}));
      connect(reaction.products[1], C.port_a)
        annotation (Line(points={{24,2},{48,2}}, color={158,66,200}));
      annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.001));
    end SimpleReaction2;

    model SimpleReaction2_Get_DfG
      "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
       extends Modelica.Icons.Example;

    Chemical.Sources.PureSubstance A
        annotation (Placement(transformation(extent={{-28,42},{-8,62}})));
      Chemical.Sensors.DissociationCoefficient reaction(nP=1, nS=2)
        annotation (Placement(transformation(extent={{10,32},{30,52}})));
    Chemical.Sources.PureSubstance B
        annotation (Placement(transformation(extent={{-28,16},{-8,36}})));

      Modelica.Blocks.Math.InverseBlockConstraints inverseBlockConstraints
        annotation (Placement(transformation(extent={{-42,-80},{82,80}})));
    Chemical.Sources.ExternalElectroChemicalPotential C(usePotentialInput=
            true)
        annotation (Placement(transformation(extent={{60,32},{40,52}})));
      Modelica.Blocks.Sources.Constant K(k=2*55.508)
        annotation (Placement(transformation(extent={{-92,-10},{-72,10}})));
    equation

      connect(reaction.DissociationCoefficient_MoleFractionBased,
        inverseBlockConstraints.u2) annotation (Line(
          points={{20,34},{20,0},{-29.6,0}},
          color={0,0,127}));
      connect(C.uInput, inverseBlockConstraints.y2) annotation (Line(
          points={{60,42},{70,42},{70,24},{46,24},{46,0},{72.7,0}},
          color={0,0,127}));
      connect(inverseBlockConstraints.u1, K.y) annotation (Line(
          points={{-48.2,0},{-71,0}},
          color={0,0,127}));
      connect(reaction.products[1], C.port_a)
        annotation (Line(points={{30,42},{40,42}}, color={158,66,200}));
      connect(A.port_a, reaction.substrates[1]) annotation (Line(points={
              {-8,52},{-8,42},{10,42},{10,41.5}}, color={158,66,200}));
      connect(B.port_a, reaction.substrates[2]) annotation (Line(points={
              {-8,26},{-8,40},{10,40},{10,42.5}}, color={158,66,200}));
      annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
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

      Chemical.Sources.PureSubstance Ag(substanceData=
            Chemical.Substances.Silver_solid())
        annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
      Chemical.Sources.PureSubstance Cl(substanceData=
            Chemical.Substances.Chloride_aqueous())
        annotation (Placement(transformation(extent={{-8,-36},{-28,-16}})));
      Chemical.Sources.PureSubstance AgCl(substanceData=
            Chemical.Substances.SilverChloride_solid())
        annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
      Chemical.Sources.ExternalIdealGasSubstance H2(
        substanceData=Chemical.Substances.Hydrogen_gas(),
        PartialPressure=100000,
        TotalPressure=100000)
        annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Chemical.Sources.PureSubstance H(substanceData=
            Chemical.Substances.Proton_aqueous())
        annotation (Placement(transformation(extent={{18,-36},{38,-16}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Chemical.Components.Reaction electrodeReaction(
        nP=2,
        p={2,2},
        nS=1) annotation (Placement(transformation(
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
          points={{-60,-18},{-42,-18},{-42,-4},{-38,-4}},
          color={158,66,200},
          thickness=1));
      connect(Cl.port_a, electrodeReaction1.substrates[2]) annotation (Line(
          points={{-28,-26},{-42,-26},{-42,-4}},
          color={158,66,200},
          thickness=1));
      connect(AgCl.port_a, electrodeReaction1.products[1]) annotation (Line(
          points={{-60,22},{-38,22},{-38,16}},
          color={158,66,200},
          thickness=1));
      connect(H.port_a, electrodeReaction.products[1]) annotation (Line(
          points={{38,-26},{50,-26},{50,-4}},
          color={158,66,200},
          thickness=1));
      connect(electrodeReaction.products[2], electrone1.port_a) annotation (Line(
          points={{54,-4},{54,-16},{66,-16}},
          color={158,66,200},
          thickness=1));
      connect(electrodeReaction1.products[2], electrone.port_a) annotation (Line(
          points={{-42,16},{-42,50},{-60,50}},
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
      connect(H2.port_a, electrodeReaction.substrates[1]) annotation (Line(
            points={{44,42},{52,42},{52,16}}, color={158,66,200}));
      annotation (
      experiment(StopTime=1), Documentation(info=
                    "<html>
<p>Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential </p>
</html>", revisions="<html>
<p><i>2015-2018</i></p>
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

      Chemical.Sources.PureSubstance Pb(substanceData=
            Chemical.Substances.Lead_solid())
        annotation (Placement(transformation(extent={{84,-34},{64,-14}})));
      Chemical.Sources.PureSubstance HSO4(substanceData=
            Chemical.Substances.HydrogenSulfate_aqueous())
        annotation (Placement(transformation(extent={{-22,-58},{-2,-38}})));
      Chemical.Sources.PureSubstance PbSO4_(substanceData=
            Chemical.Substances.LeadSulfate_solid())
        annotation (Placement(transformation(extent={{84,4},{64,24}})));
      Chemical.Sources.PureSubstance H(substanceData=
            Chemical.Substances.Proton_aqueous())
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
      Chemical.Sources.PureSubstance PbO2(substanceData=
            Chemical.Substances.LeadDioxide_solid()) annotation (Placement(
            transformation(extent={{-10,-10},{10,10}}, origin={-74,-30})));
      Chemical.Sources.PureSubstance H2O(substanceData=
            Chemical.Substances.Water_liquid())
        annotation (Placement(transformation(extent={{-2,-10},{-22,10}})));
      Chemical.Sources.PureSubstance PbSO4(substanceData=
            Chemical.Substances.LeadSulfate_solid()) annotation (Placement(
            transformation(extent={{-10,-10},{10,10}}, origin={-74,32})));

    equation
      connect(Pb.port_a, electrodeReaction1.substrates[1]) annotation (Line(
          points={{64,-24},{45.5,-24},{45.5,4},{42,4}},
          color={158,66,200},
          thickness=1));
      connect(HSO4.port_a, electrodeReaction1.substrates[2]) annotation (Line(
          points={{-2,-48},{44,-48},{44,4},{46,4}},
          color={158,66,200},
          thickness=1));
      connect(PbSO4_.port_a, electrodeReaction1.products[1]) annotation (Line(
          points={{64,14},{56,14},{56,28},{46,28},{46,24},{41.3333,24}},
          color={158,66,200},
          thickness=1));
      connect(electrodeReaction.products[1], PbSO4.port_a) annotation (Line(
          points={{-40,24},{-40,32},{-64,32}},
          color={158,66,200},
          thickness=1));
      connect(electrodeReaction.products[2], H2O.port_a) annotation (Line(
          points={{-44,24},{-40,24},{-40,32},{-34,32},{-34,0},{-22,0}},
          color={158,66,200},
          thickness=1));
      connect(PbO2.port_a, electrodeReaction.substrates[1]) annotation (Line(
          points={{-64,-30},{-42,-30},{-42,4},{-39,4}},
          color={158,66,200},
          thickness=1));
      connect(HSO4.port_a, electrodeReaction.substrates[2]) annotation (Line(
          points={{-2,-48},{-40,-48},{-40,4},{-41,4}},
          color={158,66,200},
          thickness=1));
      connect(H.port_a, electrodeReaction.substrates[3]) annotation (Line(
          points={{26,-18},{-38,-18},{-38,4},{-43,4}},
          color={158,66,200},
          thickness=1));
      connect(electrone1.port_a, electrodeReaction.substrates[4]) annotation (Line(
          points={{-66,-2},{-44,-2},{-44,4},{-45,4}},
          color={158,66,200},
          thickness=1));
      connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(
          points={{26,-18},{32,-18},{32,32},{42,32},{42,24},{44,24}},
          color={158,66,200},
          thickness=1));
      connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
          points={{64,42},{44,42},{44,24},{46.6667,24}},
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
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end StandardLeadAcidPotential;

    model StandardElectrochemicalCell2
      "Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential"
     extends Modelica.Icons.Example;

      Chemical.Components.Solution cathode(ElectricGround=false)
        annotation (Placement(transformation(extent={{-90,-40},{-46,68}})));

      Chemical.Components.Solution anode(ElectricGround=false)
        annotation (Placement(transformation(extent={{60,-40},{96,70}})));

      Chemical.Sources.PureSubstance H2O(substanceData=
            Chemical.Substances.Water_liquid())
        annotation (Placement(transformation(extent={{-8,-36},{-28,-16}})));
      Chemical.Sources.PureSubstance O2(redeclare package stateOfMatter =
            Interfaces.IdealGas, substanceData=Chemical.Substances.Oxygen_gas())
        annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
      Chemical.Sources.ExternalIdealGasSubstance H2(
        substanceData=Chemical.Substances.Hydrogen_gas(),
        PartialPressure=100000,
        TotalPressure=100000)
        annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Chemical.Sources.PureSubstance H(substanceData=
            Chemical.Substances.Proton_aqueous())
        annotation (Placement(transformation(extent={{18,-36},{38,-16}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Chemical.Components.Reaction electrodeReaction(
        nP=2,
        p={2,2},
        nS=1) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={52,6})));
      Chemical.Components.Reaction electrodeReaction1(
        s={4},
        p={2,8,8},                                    nS=1,
        nP=3)                                                     annotation (
         Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-40,6})));
    Chemical.Components.ElectronTransfer electrone
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    Chemical.Components.ElectronTransfer electrone1
        annotation (Placement(transformation(extent={{86,-26},{66,-6}})));

    equation
      connect(H2O.port_a, electrodeReaction1.substrates[1]) annotation (Line(
          points={{-28,-26},{-40,-26},{-40,-4}},
          color={158,66,200},
          thickness=1));
      connect(H.port_a, electrodeReaction.products[1]) annotation (Line(
          points={{38,-26},{50,-26},{50,-4}},
          color={158,66,200},
          thickness=1));
      connect(electrodeReaction.products[2], electrone1.port_a) annotation (Line(
          points={{54,-4},{54,-16},{66,-16}},
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
      connect(H2.port_a, electrodeReaction.substrates[1]) annotation (Line(
            points={{44,42},{52,42},{52,16}}, color={158,66,200}));
      connect(O2.port_a, electrodeReaction1.products[1]) annotation (Line(
            points={{-60,22},{-34,22},{-34,16},{-37.3333,16}}, color={158,66,
              200}));
      connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(points=
             {{38,-26},{12,-26},{12,30},{-40,30},{-40,16}}, color={158,66,200}));
      connect(electrone.port_a, electrodeReaction1.products[3]) annotation (
          Line(points={{-60,50},{-42,50},{-42,16},{-42.6667,16}}, color={158,66,
              200}));
      annotation (
      experiment(StopTime=1), Documentation(info=
                    "<html>
<p>Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential </p>
</html>", revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end StandardElectrochemicalCell2;
  end CheckSubstancesData;

  model WaterElectrolysis "Water electrolysis"
    extends Modelica.Icons.Example;
    Chemical.Components.Substance O2_gas(
      substanceData=Chemical.Substances.Oxygen_gas(),
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      use_mass_start=false,
      amountOfSubstance_start=0.001)
      annotation (Placement(transformation(extent={{-22,-6},{-2,14}})));

    Chemical.Components.Substance H2_gas(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Hydrogen_gas(),
      use_mass_start=false,
      amountOfSubstance_start=0.001)
      annotation (Placement(transformation(extent={{16,-6},{36,14}})));
    Chemical.Components.Solution anode(ElectricGround=false)
      annotation (Placement(transformation(extent={{58,-78},{92,30}})));
    Chemical.Components.Solution cathode(ElectricGround=false)
      annotation (Placement(transformation(extent={{-90,-80},{-56,28}})));
    Chemical.Components.Solution water(temperature_start=310.15)
      annotation (Placement(transformation(extent={{-28,-80},{18,-46}})));
    Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
      annotation (Placement(transformation(extent={{-42,70},{-22,90}})));
    Chemical.Components.Reaction reaction(
      nS=2,
      s={2,4},
      nP=3,
      p={2,1,4}) annotation (Placement(transformation(
          extent={{11,11},{-11,-11}},
          rotation=180,
          origin={-23,-31})));
  Chemical.Components.ElectronTransfer electrone
      annotation (Placement(transformation(extent={{84,-38},{64,-18}})));
  Chemical.Components.ElectronTransfer electrone1
      annotation (Placement(transformation(extent={{-84,-34},{-64,-14}})));
  Modelica.Electrical.Analog.Basic.Resistor resistor(R=1)
    annotation (Placement(transformation(extent={{-36,38},{-16,58}})));
  Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
    annotation (Placement(transformation(extent={{-66,38},{-46,58}})));
    Chemical.Components.Solution air(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGas)
      annotation (Placement(transformation(extent={{-40,-16},{50,26}})));
    Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
      annotation (Placement(transformation(extent={{18,38},{-2,58}})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{38,26},{58,46}})));
    Components.Substance liquidWater(substanceData=
          Chemical.Substances.Water_liquid(), mass_start=1)
      annotation (Placement(transformation(extent={{-14,-72},{6,-52}})));
  equation
  connect(electrone1.pin,voltageSensor. p) annotation (Line(
      points={{-84,-24},{-92,-24},{-92,48},{-74,48},{-74,80},{-42,80}},
      color={0,0,255}));
  connect(electrone.pin,voltageSensor. n) annotation (Line(
      points={{84,-28},{84,48},{60,48},{60,80},{-22,80}},
      color={0,0,255}));
  connect(electrone.solution,anode. solution) annotation (Line(
      points={{80,-38},{80,-76.92},{85.2,-76.92}},
      color={127,127,0}));
  connect(electrone1.pin,currentSensor. p) annotation (Line(
      points={{-84,-24},{-92,-24},{-92,48},{-66,48}},
      color={0,0,255}));
  connect(currentSensor.n,resistor. p) annotation (Line(
      points={{-46,48},{-36,48}},
      color={0,0,255}));
  connect(electrone1.solution,cathode. solution) annotation (Line(
      points={{-80,-34},{-80,-66},{-74,-66},{-74,-78.92},{-62.8,-78.92}},
      color={127,127,0}));
    connect(O2_gas.solution, air.solution) annotation (Line(points={{-18,-6},{-16,-6},
            {-16,-15.58},{32,-15.58}},     color={127,127,0}));
    connect(H2_gas.solution, air.solution) annotation (Line(points={{20,-6},{20,-15.58},
            {32,-15.58}}, color={127,127,0}));
    connect(electrone1.port_a, reaction.substrates[2]) annotation (Line(points={{-64,-24},
            {-48,-24},{-48,-33.2},{-34,-33.2}},      color={158,66,200}));
    connect(H2_gas.port_a, reaction.products[1]) annotation (Line(points={{36,4},{42,4},
            {42,-34},{-12,-34},{-12,-28.0667}},       color={158,66,200}));
    connect(O2_gas.port_a, reaction.products[2]) annotation (Line(points={{-2,4},{
            2,4},{2,-31},{-12,-31}}, color={158,66,200}));
    connect(electrone.port_a, reaction.products[3]) annotation (Line(points={{64,-28},
            {26,-28},{26,-33.9333},{-12,-33.9333}}, color={158,66,200}));
    connect(constantVoltage.p, voltageSensor.n) annotation (Line(points={{18,48},{
            60,48},{60,80},{-22,80}}, color={0,0,255}));
    connect(resistor.n, constantVoltage.n)
      annotation (Line(points={{-16,48},{-2,48}}, color={0,0,255}));
    connect(constantVoltage.p, ground.p)
      annotation (Line(points={{18,48},{48,48},{48,46}}, color={0,0,255}));
    connect(liquidWater.solution, water.solution) annotation (Line(points={{-10,-72},{
            -10,-79.66},{8.8,-79.66}},       color={127,127,0}));
    connect(liquidWater.port_a, reaction.substrates[1]) annotation (Line(points={{6,-62},
            {-48,-62},{-48,-28.8},{-34,-28.8}},         color={158,66,200}));
    annotation ( experiment(StopTime=1), Documentation(info="<html>
<p>The water ecectrolysis: </p>
<p><b>2 H<sub>2</sub>O +&nbsp;&nbsp;4 e<sup>-</sup><sub>(catode)</sub>&nbsp;&lt;-&gt;  2 H<sub>2</sub> + O<sub>2</sub>&nbsp;+&nbsp;&nbsp;4 e<sup>-</sup><sub>(anode)</sub>&nbsp;</b></p>
</html>"));
  end WaterElectrolysis;

  model H2O_ElectrochemicalCell
   extends Modelica.Icons.Example;

    Chemical.Components.Solution cathode(ElectricGround=false)
      annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
    Chemical.Components.Solution anode(ElectricGround=false)
      annotation (Placement(transformation(extent={{62,-50},{96,50}})));

    Chemical.Components.Solution solution1(ElectricGround=false)
      annotation (Placement(transformation(extent={{-30,-60},{38,6}})));

    Chemical.Sources.ExternalIdealGasSubstance H2(
      substanceData=Chemical.Substances.Hydrogen_gas(),
      PartialPressure=100000,
      TotalPressure=100000)
      annotation (Placement(transformation(extent={{24,32},{44,52}})));
    Chemical.Components.Substance H(
      substanceData=Chemical.Substances.Proton_aqueous(),
      use_mass_start=false,
      amountOfSubstance_start=1e-7)
      annotation (Placement(transformation(extent={{-6,-22},{14,-2}})));
    Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
      annotation (Placement(transformation(extent={{-6,64},{14,84}})));
    Chemical.Components.Reaction electrodeReaction(
      p={2,2},
      nS=1,
      nP=2) annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={52,6})));
    Chemical.Components.Reaction electrodeReaction1(
      s={4},
      p={2,8,8},
      nS=1,
      nP=3)                                                     annotation (
        Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=90,
          origin={-40,0})));

  Chemical.Components.ElectronTransfer electrone
      annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                               //(substanceData=Chemical.Examples.Substances.Electrone_solid())
  Chemical.Components.ElectronTransfer electrone1
      annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                //(substanceData=Chemical.Examples.Substances.Electrone_solid())
  Modelica.Electrical.Analog.Basic.Ground ground
    annotation (Placement(transformation(extent={{62,54},{82,74}})));
    Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
        mass_start=1)
      annotation (Placement(transformation(extent={{-6,-54},{14,-34}})));
    Sources.ExternalIdealGasSubstance O2_(
      substanceData=Chemical.Substances.Oxygen_gas(),
      PartialPressure=100000,
      TotalPressure=100000)
      annotation (Placement(transformation(extent={{-6,32},{-26,52}})));
  equation
    connect(H.solution, solution1.solution) annotation (Line(points={{-2,-22},{
            -2,-30},{24.4,-30},{24.4,-59.34}},
                                       color={127,127,0}));
  connect(electrone.solution, cathode.solution) annotation (Line(
      points={{-74,32},{-74,-18},{-64,-18},{-64,-42.84},{-54.4,-42.84}},
      color={127,127,0}));
  connect(electrone1.solution, anode.solution) annotation (Line(
      points={{84,-26},{84,-49},{89.2,-49}},
      color={127,127,0}));
    connect(voltageSensor.p, electrone.pin) annotation (Line(
        points={{-6,74},{-96,74},{-96,42},{-78,42}},
        color={0,0,255}));
    connect(voltageSensor.n, electrone1.pin) annotation (Line(
        points={{14,74},{92,74},{92,-16},{88,-16}},
        color={0,0,255}));
    connect(electrone1.pin, ground.p) annotation (Line(
        points={{88,-16},{92,-16},{92,74},{72,74}},
        color={0,0,255}));
    connect(H2O.solution, solution1.solution) annotation (Line(points={{-2,-54},
            {-2,-59.34},{24.4,-59.34}}, color={127,127,0}));
    connect(H2.port_a, electrodeReaction.substrates[1])
      annotation (Line(points={{44,42},{52,42},{52,16}}, color={158,66,200}));
    connect(electrodeReaction1.substrates[1], H2O.port_a) annotation (Line(
          points={{-40,-10},{-40,-44},{14,-44}}, color={158,66,200}));
    connect(O2_.port_a, electrodeReaction1.products[1]) annotation (Line(points={{-26,42},
            {-38,42},{-38,10},{-37.3333,10}},          color={158,66,200}));
    connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(points={
            {14,-12},{20,-12},{20,20},{-40,20},{-40,10}}, color={158,66,200}));
    connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
          points={{-58,42},{-42,42},{-42,10},{-42.6667,10}}, color={158,66,200}));
    connect(H.port_a, electrodeReaction.products[1]) annotation (Line(points={{
            14,-12},{50,-12},{50,-4}}, color={158,66,200}));
    connect(electrone1.port_a, electrodeReaction.products[2]) annotation (Line(
          points={{68,-16},{54,-16},{54,-4}}, color={158,66,200}));
    annotation (
    experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end H2O_ElectrochemicalCell;

  package ClimateChange
    class References "References"
      extends Modelica.Icons.References;

      annotation (DocumentationClass=true,Documentation(info="<html>
<table cellspacing=\"0\" cellpadding=\"2\" border=\"0\">
<tr>
<td>[Nelabhotla2019]</td>
<td>Anirudh Bhanu Teja Nelabhotla, Rune Bakke, Carlos Dinamarca,
	\"Performance Analysis of Biocathode in Bioelectrochemical CO2 Reduction\"
	 Catalysts, 9, 683, 2019,
	<a href=\"https://doi.org/10.3390/catal9080683\">doi:10.3390/catal9080683</a>.
</td>
</tr>
</table>
</html>"));
    end References;

    model AcetoclasticMethanogenesis
      extends Modelica.Icons.Example;
      Components.Solution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
      Components.Substance CH3COOH(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.AceticAcid_aqueous(),
        mass_start=0.001) "Acetic acid"
        annotation (Placement(transformation(extent={{-72,30},{-52,50}})));
      Components.Substance CH4(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.Methan_aqueous(),
        mass_start=0.001) "Methan"
        annotation (Placement(transformation(extent={{26,48},{46,68}})));
      Components.Substance CO2(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
        mass_start=0.001)
        annotation (Placement(transformation(extent={{56,10},{76,30}})));
      Components.Substance Water(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.Water_liquid(),
        mass_start=1)
        annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));
      Components.Reaction reaction(nS=1, nP=2)
        "Acetoclastic (heterotrophic) methanogenesis"
        annotation (Placement(transformation(extent={{-20,30},{0,50}})));
    equation
      connect(CH3COOH.port_a, reaction.substrates[1])
        annotation (Line(points={{-52,40},{-20,40}}, color={158,66,200}));
      connect(CH4.port_a, reaction.products[1]) annotation (Line(points={{46,58},
              {24,58},{24,42},{0,42}}, color={158,66,200}));
      connect(CO2.port_a, reaction.products[2]) annotation (Line(points={{76,20},
              {24,20},{24,38},{0,38}}, color={158,66,200}));
      connect(CH3COOH.solution, solution.solution) annotation (Line(points={{
              -68,30},{-68,-88},{60,-88},{60,-98}}, color={127,127,0}));
      connect(Water.solution, solution.solution) annotation (Line(points={{-10,
              -66},{-10,-88},{60,-88},{60,-98}}, color={127,127,0}));
      connect(CO2.solution, solution.solution)
        annotation (Line(points={{60,10},{60,-98}}, color={127,127,0}));
      connect(CH4.solution, solution.solution) annotation (Line(points={{30,48},
              {32,48},{32,-88},{60,-88},{60,-98}}, color={127,127,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end AcetoclasticMethanogenesis;

    model HydrogenotrophicMethanogenesis
      extends Modelica.Icons.Example;
      Components.Solution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
      Components.Substance CH4(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.Methan_aqueous(),
        mass_start=0.001) "Methan"
        annotation (Placement(transformation(extent={{32,70},{52,90}})));
      Components.Substance CO2(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
        mass_start=0.001)
        annotation (Placement(transformation(extent={{-68,40},{-48,60}})));
      Components.Substance H2O(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.Water_liquid(),
        mass_start=1)
        annotation (Placement(transformation(extent={{36,36},{56,56}})));
      Components.Reaction reaction(
        KC=1e-7,
        s={4,1},
        p={1,2},
        nS=2,
        nP=2) "Hydrogenotrophic (autotrophic) methanogenesis"
        annotation (Placement(transformation(extent={{-18,50},{2,70}})));
      Components.Substance H2(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.Hydrogen_aqueous(),
        mass_start=0.001)
        annotation (Placement(transformation(extent={{-74,66},{-54,86}})));
    equation
      connect(H2O.solution, solution.solution) annotation (Line(points={{40,36},
              {40,22},{60,22},{60,-98}}, color={127,127,0}));
      connect(CO2.solution, solution.solution) annotation (Line(points={{-64,40},
              {-64,22},{60,22},{60,-98}}, color={127,127,0}));
      connect(CH4.solution, solution.solution) annotation (Line(points={{36,70},
              {36,22},{60,22},{60,-98}}, color={127,127,0}));
      connect(H2.solution, solution.solution) annotation (Line(points={{-70,66},
              {-70,22},{60,22},{60,-98}}, color={127,127,0}));
      connect(H2.port_a, reaction.substrates[1]) annotation (Line(points={{-54,
              76},{-36,76},{-36,62},{-18,62}}, color={158,66,200}));
      connect(CO2.port_a, reaction.substrates[2]) annotation (Line(points={{-48,
              50},{-32,50},{-32,58},{-18,58}}, color={158,66,200}));
      connect(CH4.port_a, reaction.products[1]) annotation (Line(points={{52,80},
              {26,80},{26,62},{2,62}}, color={158,66,200}));
      connect(H2O.port_a, reaction.products[2]) annotation (Line(points={{56,46},
              {30,46},{30,58},{2,58}}, color={158,66,200}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=1,
          Tolerance=1e-08,
          __Dymola_Algorithm="Dassl"));
    end HydrogenotrophicMethanogenesis;

    model MethanElectrosynthesis
      "Direct electron transfer (electrosynthesis reaction-bioelectrochemical methane)"
     extends Modelica.Icons.Example;

      Chemical.Components.Solution cathode(ElectricGround=false)
        annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
      Chemical.Components.Solution anode(ElectricGround=false)
        annotation (Placement(transformation(extent={{62,-50},{96,50}})));

      Chemical.Components.Solution solution1(ElectricGround=false)
        annotation (Placement(transformation(extent={{-30,-96},{38,6}})));

      Chemical.Sources.ExternalIdealGasSubstance CH4(
        substanceData=Chemical.Substances.Methan_gas(),
        PartialPressure=100000,
        TotalPressure=100000)
        annotation (Placement(transformation(extent={{22,26},{42,46}})));
      Chemical.Components.Substance H(
        substanceData=Chemical.Substances.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7)
        annotation (Placement(transformation(extent={{-6,-22},{14,-2}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-20,62},{0,82}})));
      Chemical.Components.Reaction electrodeReaction(
        s={1,8,8},
        p={1,2},
        nP=2,
        nS=3) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={52,6})));
      Chemical.Components.Reaction electrodeReaction1(
        s={4},
        p={2,8,8},
        nS=1,
        nP=3)                                                     annotation (
          Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-40,0})));

    Chemical.Components.ElectronTransfer electrone
        annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                                 //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Chemical.Components.ElectronTransfer electrone1
        annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                  //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{46,52},{66,72}})));
      Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
          mass_start=1)
        annotation (Placement(transformation(extent={{-6,-80},{14,-60}})));
      Sources.ExternalIdealGasSubstance O2_(
        substanceData=Chemical.Substances.Oxygen_gas(),
        PartialPressure=100000,
        TotalPressure=100000)
        annotation (Placement(transformation(extent={{-6,32},{-26,52}})));
      Sources.ExternalIdealGasSubstance CO2(
        substanceData=Chemical.Substances.CarbonDioxide_gas(),
        PartialPressure=100000,
        TotalPressure=100000)
        annotation (Placement(transformation(extent={{88,-90},{68,-70}})));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
        annotation (Placement(transformation(extent={{-80,78},{-60,98}})));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R=0.01)
        annotation (Placement(transformation(extent={{12,78},{32,98}})));
      Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor
        annotation (Placement(transformation(extent={{76,78},{96,98}})));
    equation
      connect(H.solution, solution1.solution) annotation (Line(points={{-2,-22},
              {-2,-30},{24.4,-30},{24.4,-94.98}},
                                         color={127,127,0}));
    connect(electrone.solution, cathode.solution) annotation (Line(
        points={{-74,32},{-74,-18},{-64,-18},{-64,-42.84},{-54.4,-42.84}},
        color={127,127,0}));
    connect(electrone1.solution, anode.solution) annotation (Line(
        points={{84,-26},{84,-49},{89.2,-49}},
        color={127,127,0}));
      connect(voltageSensor.p, electrone.pin) annotation (Line(
          points={{-20,72},{-88,72},{-88,42},{-78,42}},
          color={0,0,255}));
      connect(voltageSensor.n, electrone1.pin) annotation (Line(
          points={{0,72},{98,72},{98,-16},{88,-16}},
          color={0,0,255}));
      connect(electrone1.pin, ground.p) annotation (Line(
          points={{88,-16},{98,-16},{98,72},{56,72}},
          color={0,0,255}));
      connect(H2O.solution, solution1.solution) annotation (Line(points={{-2,-80},
              {-2,-94.98},{24.4,-94.98}}, color={127,127,0}));
      connect(electrodeReaction1.substrates[1], H2O.port_a) annotation (Line(
            points={{-40,-10},{-40,-70},{14,-70}}, color={158,66,200}));
      connect(O2_.port_a, electrodeReaction1.products[1]) annotation (Line(
            points={{-26,42},{-38,42},{-38,10},{-37.3333,10}}, color={158,66,
              200}));
      connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(points=
             {{14,-12},{20,-12},{20,20},{-40,20},{-40,10}}, color={158,66,200}));
      connect(electrone.port_a, electrodeReaction1.products[3]) annotation (
          Line(points={{-58,42},{-42,42},{-42,10},{-42.6667,10}}, color={158,66,
              200}));
      connect(CH4.port_a, electrodeReaction.products[1]) annotation (Line(
            points={{42,36},{54,36},{54,16}}, color={158,66,200}));
      connect(H2O.port_a, electrodeReaction.products[2]) annotation (Line(
            points={{14,-70},{42,-70},{42,24},{50,24},{50,16}}, color={158,66,
              200}));
      connect(CO2.port_a, electrodeReaction.substrates[1]) annotation (Line(
            points={{68,-80},{54.6667,-80},{54.6667,-4}}, color={158,66,200}));
      connect(H.port_a, electrodeReaction.substrates[2]) annotation (Line(
            points={{14,-12},{52,-12},{52,-4}}, color={158,66,200}));
      connect(electrone1.port_a, electrodeReaction.substrates[3]) annotation (
          Line(points={{68,-16},{50,-16},{50,-4},{49.3333,-4}}, color={158,66,
              200}));
      connect(electrone.pin, constantVoltage.p) annotation (Line(points={{-78,42},
              {-88,42},{-88,88},{-80,88}},       color={0,0,255}));
      connect(constantVoltage.n, resistor.p)
        annotation (Line(points={{-60,88},{12,88}},   color={0,0,255}));
      connect(powerSensor.nv, resistor.n) annotation (Line(points={{86,78},{74,
              78},{74,88},{32,88}},   color={0,0,255}));
      connect(powerSensor.pv, constantVoltage.p) annotation (Line(points={{86,98},
              {-82,98},{-82,88},{-80,88}},         color={0,0,255}));
      connect(resistor.n, powerSensor.pc) annotation (Line(points={{32,88},{76,
              88}},                     color={0,0,255}));
      connect(powerSensor.nc, electrone1.pin) annotation (Line(points={{96,88},
              {98,88},{98,-16},{88,-16}},    color={0,0,255}));
      annotation (
      experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MethanElectrosynthesis;

    model ElectrolysisMethanation
     extends Modelica.Icons.Example;

      Chemical.Components.Solution cathode(ElectricGround=false)
        annotation (Placement(transformation(extent={{-86,26},{-46,72}})));
      Chemical.Components.Solution anode(ElectricGround=false)
        annotation (Placement(transformation(extent={{58,26},{94,72}})));

      Chemical.Components.Solution solution1(ElectricGround=false)
        annotation (Placement(transformation(extent={{-90,-90},{92,14}})));

      Chemical.Components.Substance H(
        substanceData=Chemical.Substances.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7)
        annotation (Placement(transformation(extent={{-84,-14},{-64,6}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Chemical.Components.Reaction electrodeReaction(
        p={2,2},
        nP=2,
        nS=1) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={50,16})));
      Chemical.Components.Reaction electrodeReaction1(
        s={4},
        p={2,8,8},
        nS=1,
        nP=3)                                                     annotation (
          Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-36,28})));

    Chemical.Components.ElectronTransfer electrone
        annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
                                 //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Chemical.Components.ElectronTransfer electrone1
        annotation (Placement(transformation(extent={{88,38},{68,58}})));
                                  //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{22,54},{42,74}})));
      Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
          mass_start=1)
        annotation (Placement(transformation(extent={{78,-28},{58,-8}})));
      Sources.ExternalIdealGasSubstance O2_(
        substanceData=Chemical.Substances.Oxygen_gas(),
        PartialPressure=100000,
        TotalPressure=100000)
        annotation (Placement(transformation(extent={{0,36},{-20,56}})));
      Components.Substance CH4(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.Methan_aqueous(),
        mass_start=0.001) "Methan"
        annotation (Placement(transformation(extent={{74,-64},{54,-44}})));
      Components.Substance CO2(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
        mass_start=0.1)
        annotation (Placement(transformation(extent={{-72,-76},{-52,-56}})));
      Components.Reaction reaction(
        KC=1e-7,
        s={4,1},
        p={1,2},
        nS=2,
        nP=2) "Hydrogenotrophic (autotrophic) methanogenesis"
        annotation (Placement(transformation(extent={{0,-54},{20,-34}})));
      Components.Substance H2(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.Hydrogen_aqueous(),
        mass_start=5e-11)
        annotation (Placement(transformation(extent={{-76,-44},{-56,-24}})));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
        annotation (Placement(transformation(extent={{-72,84},{-52,104}})));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R=0.01)
        annotation (Placement(transformation(extent={{20,84},{40,104}})));
      Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor
        annotation (Placement(transformation(extent={{112,80},{132,100}})));
    equation
      connect(H.solution, solution1.solution) annotation (Line(points={{-80,-14},{-80,
              -88},{55.6,-88},{55.6,-88.96}},
                                         color={127,127,0}));
    connect(electrone.solution, cathode.solution) annotation (Line(
        points={{-74,44},{-54,44},{-54,26.46}},
        color={127,127,0}));
    connect(electrone1.solution, anode.solution) annotation (Line(
        points={{84,38},{84,26.46},{86.8,26.46}},
        color={127,127,0}));
      connect(voltageSensor.p, electrone.pin) annotation (Line(
          points={{-6,74},{-96,74},{-96,54},{-78,54}},
          color={0,0,255}));
      connect(voltageSensor.n, electrone1.pin) annotation (Line(
          points={{14,74},{92,74},{92,48},{88,48}},
          color={0,0,255}));
      connect(electrone1.pin, ground.p) annotation (Line(
          points={{88,48},{92,48},{92,74},{32,74}},
          color={0,0,255}));
      connect(H2O.solution, solution1.solution) annotation (Line(points={{74,-28},{74,
              -88.96},{55.6,-88.96}},     color={127,127,0}));
      connect(electrodeReaction1.substrates[1], H2O.port_a) annotation (Line(
            points={{-36,18},{-36,4},{36,4},{36,-18},{58,-18}},
                                                   color={158,66,200}));
      connect(O2_.port_a, electrodeReaction1.products[1]) annotation (Line(points={{-20,46},
              {-38,46},{-38,38},{-33.3333,38}},         color={158,66,200}));
      connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(points={{-64,
              -4},{18,-4},{18,62},{-36,62},{-36,38}}, color={158,66,200}));
      connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
            points={{-58,54},{-42,54},{-42,38},{-38.6667,38}}, color={158,66,200}));
      connect(H.port_a, electrodeReaction.products[1])
        annotation (Line(points={{-64,-4},{48,-4},{48,6}}, color={158,66,200}));
      connect(electrone1.port_a, electrodeReaction.products[2]) annotation (Line(
            points={{68,48},{60,48},{60,-2},{52,-2},{52,6}}, color={158,66,200}));
      connect(H2.port_a, reaction.substrates[1]) annotation (Line(points={{-56,-34},
              {-30,-34},{-30,-42},{0,-42}}, color={158,66,200}));
      connect(CO2.port_a, reaction.substrates[2]) annotation (Line(points={{-52,-66},
              {-30,-66},{-30,-46},{0,-46}}, color={158,66,200}));
      connect(CH4.port_a, reaction.products[1]) annotation (Line(points={{54,-54},
              {40,-54},{40,-46},{20,-46},{20,-42}},  color={158,66,200}));
      connect(H2.port_a, electrodeReaction.substrates[1]) annotation (Line(points={{
              -56,-34},{-30,-34},{-30,32},{50,32},{50,26}}, color={158,66,200}));
      connect(reaction.products[2], H2O.port_a) annotation (Line(points={{20,-46},
              {20,-42},{40,-42},{40,-18},{58,-18}},
                                                color={158,66,200}));
      connect(CO2.solution, solution1.solution) annotation (Line(points={{-68,-76},{
              -80,-76},{-80,-88.96},{55.6,-88.96}}, color={127,127,0}));
      connect(H2.solution, solution1.solution) annotation (Line(points={{-72,-44},{-80,
              -44},{-80,-88},{56,-88},{56,-88.96},{55.6,-88.96}}, color={127,127,0}));
      connect(CH4.solution, solution1.solution) annotation (Line(points={{70,-64},{74,
              -64},{74,-88.96},{55.6,-88.96}}, color={127,127,0}));
      connect(electrone.pin, constantVoltage.p) annotation (Line(points={{-78,
              54},{-96,54},{-96,94},{-72,94}}, color={0,0,255}));
      connect(constantVoltage.n, resistor.p)
        annotation (Line(points={{-52,94},{20,94}}, color={0,0,255}));
      connect(powerSensor.nv, resistor.n) annotation (Line(points={{122,80},{82,
              80},{82,94},{40,94}}, color={0,0,255}));
      connect(powerSensor.pv, constantVoltage.p) annotation (Line(points={{122,
              100},{-74,100},{-74,94},{-72,94}}, color={0,0,255}));
      connect(resistor.n, powerSensor.pc) annotation (Line(points={{40,94},{76,
              94},{76,90},{112,90}}, color={0,0,255}));
      connect(powerSensor.nc, electrone1.pin) annotation (Line(points={{132,90},
              {142,90},{142,48},{88,48}}, color={0,0,255}));
      annotation (
      experiment(StopTime=0.0001, __Dymola_Algorithm="Dassl"),
                                   Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ElectrolysisMethanation;

    model ElectrochemicalAcetateProduction
      "Direct electron transfer (electrochemical acetate production)"
     extends Modelica.Icons.Example;

      Chemical.Components.Solution cathode(ElectricGround=false)
        annotation (Placement(transformation(extent={{-86,26},{-46,72}})));
      Chemical.Components.Solution anode(ElectricGround=false)
        annotation (Placement(transformation(extent={{58,26},{94,72}})));

      Chemical.Components.Solution solution1(ElectricGround=false)
        annotation (Placement(transformation(extent={{-90,-90},{92,14}})));

      Chemical.Components.Substance H(
        substanceData=Chemical.Substances.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-4)
        annotation (Placement(transformation(extent={{-84,-14},{-64,6}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Chemical.Components.Reaction electrodeReaction1(
        s={4},
        p={2,8,8},
        nS=1,
        nP=3)                                                     annotation (
          Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-36,28})));

    Chemical.Components.ElectronTransfer electrone
        annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
                                 //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Chemical.Components.ElectronTransfer electrone1
        annotation (Placement(transformation(extent={{88,38},{68,58}})));
                                  //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{22,54},{42,74}})));
      Components.Substance H2O(substanceData=Chemical.Substances.Water_liquid(),
          mass_start=1)
        annotation (Placement(transformation(extent={{78,-28},{58,-8}})));
      Sources.ExternalIdealGasSubstance O2_(
        substanceData=Chemical.Substances.Oxygen_gas(),
        PartialPressure=100000,
        TotalPressure=100000)
        annotation (Placement(transformation(extent={{0,36},{-20,56}})));
      Components.Substance AcAc(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.AceticAcid_aqueous(),
        mass_start=0.001) "Acetic Acid"
        annotation (Placement(transformation(extent={{60,-60},{40,-40}})));
      Components.Substance CO2(
        redeclare package stateOfMatter = Interfaces.Incompressible,
        substanceData=Chemical.Substances.CarbonDioxide_aqueous(),
        mass_start=1)
        annotation (Placement(transformation(extent={{-10,-34},{10,-14}})));
      Components.Reaction reaction(
        s={2,8,8},
        p={1,2},
        nS=3,
        nP=2) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={50,30})));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
        annotation (Placement(transformation(extent={{-54,76},{-34,96}})));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R=0.01)
        annotation (Placement(transformation(extent={{38,76},{58,96}})));
      Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor
        annotation (Placement(transformation(extent={{130,72},{150,92}})));
    equation
      connect(H.solution, solution1.solution) annotation (Line(points={{-80,-14},{-80,
              -88},{55.6,-88},{55.6,-88.96}},
                                         color={127,127,0}));
    connect(electrone.solution, cathode.solution) annotation (Line(
        points={{-74,44},{-54,44},{-54,26.46}},
        color={127,127,0}));
    connect(electrone1.solution, anode.solution) annotation (Line(
        points={{84,38},{84,26.46},{86.8,26.46}},
        color={127,127,0}));
      connect(voltageSensor.p, electrone.pin) annotation (Line(
          points={{-6,74},{-96,74},{-96,54},{-78,54}},
          color={0,0,255}));
      connect(voltageSensor.n, electrone1.pin) annotation (Line(
          points={{14,74},{92,74},{92,48},{88,48}},
          color={0,0,255}));
      connect(electrone1.pin, ground.p) annotation (Line(
          points={{88,48},{92,48},{92,74},{32,74}},
          color={0,0,255}));
      connect(H2O.solution, solution1.solution) annotation (Line(points={{74,-28},{74,
              -88.96},{55.6,-88.96}},     color={127,127,0}));
      connect(electrodeReaction1.substrates[1], H2O.port_a) annotation (Line(
            points={{-36,18},{-36,2},{38,2},{38,-18},{58,-18}},
                                                   color={158,66,200}));
      connect(O2_.port_a, electrodeReaction1.products[1]) annotation (Line(points={{-20,46},
              {-38,46},{-38,38},{-33.3333,38}},         color={158,66,200}));
      connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(points={{-64,-4},
              {-28,-4},{-28,58},{-36,58},{-36,38}},   color={158,66,200}));
      connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
            points={{-58,54},{-42,54},{-42,38},{-38.6667,38}}, color={158,66,200}));
      connect(CO2.solution, solution1.solution) annotation (Line(points={{-6,-34},
              {-80,-34},{-80,-88.96},{55.6,-88.96}},color={127,127,0}));
      connect(AcAc.solution, solution1.solution) annotation (Line(points={{56,
              -60},{58,-60},{58,-88.96},{55.6,-88.96}}, color={127,127,0}));
      connect(CO2.port_a, reaction.substrates[1]) annotation (Line(points={{10,-24},
              {20,-24},{20,52},{52.6667,52},{52.6667,40}},      color={158,66,
              200}));
      connect(H.port_a, reaction.substrates[2]) annotation (Line(points={{-64,
              -4},{14,-4},{14,56},{50,56},{50,40}}, color={158,66,200}));
      connect(electrone1.port_a, reaction.substrates[3]) annotation (Line(
            points={{68,48},{46,48},{46,40},{47.3333,40}}, color={158,66,200}));
      connect(AcAc.port_a, reaction.products[1]) annotation (Line(points={{40,
              -50},{24,-50},{24,8},{52,8},{52,20}}, color={158,66,200}));
      connect(H2O.port_a, reaction.products[2]) annotation (Line(points={{58,
              -18},{58,-19},{48,-19},{48,20}}, color={158,66,200}));
      connect(electrone.pin, constantVoltage.p) annotation (Line(points={{-78,
              54},{-96,54},{-96,86},{-54,86}}, color={0,0,255}));
      connect(constantVoltage.n, resistor.p)
        annotation (Line(points={{-34,86},{38,86}}, color={0,0,255}));
      connect(powerSensor.nv, resistor.n) annotation (Line(points={{140,72},{
              100,72},{100,86},{58,86}}, color={0,0,255}));
      connect(powerSensor.pv, constantVoltage.p) annotation (Line(points={{140,
              92},{-56,92},{-56,86},{-54,86}}, color={0,0,255}));
      connect(resistor.n, powerSensor.pc) annotation (Line(points={{58,86},{94,
              86},{94,82},{130,82}}, color={0,0,255}));
      connect(powerSensor.nc, electrone1.pin) annotation (Line(points={{150,82},
              {162,82},{162,48},{88,48}}, color={0,0,255}));
      annotation (
      experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ElectrochemicalAcetateProduction;
  end ClimateChange;
end Examples;
