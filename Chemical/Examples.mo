within Chemical;
package Examples "Tests for top level components of undirected"
  extends Modelica.Icons.ExamplesPackage;

  model Definitions
    extends Modelica.Icons.Example;

    import Chemical.Interfaces.Definition;
    import Chemical.Substances.Gas;
    import Chemical.Substances.Liquid;
    import Chemical.Interfaces.processData;
    constant Real R = Modelica.Constants.R;

    constant Definition O2_aq =
      Gas.O2 + processData(0.0013,-1500*R);

    constant Definition H2O_formation =
      Gas.H2O - (Gas.H2 + 0.5*Gas.O2);

    constant Definition Hemoglobin =
      Liquid.Unknown;

    import Chemical.Interfaces.SolutionState;
    import Chemical.Interfaces.Phase;

    constant SolutionState SATP =
      SolutionState(phase=Phase.Gas,
                  T=298.15);

    SolutionState heatingSolution =
      SolutionState(phase=Phase.Gas,
                  T=273.15+time);

    import Chemical.Interfaces.ProcessProperties;

    ProcessProperties O2_dissolving_props(
       definition=processData(0.0013,-1500*R),
       solutionState=heatingSolution);

    ProcessProperties H2O_formation_props(
       definition=H2O_formation,
       solutionState=heatingSolution);

    ProcessProperties H2O_vaporization_props(
       definition=Gas.H2O - Liquid.H2O,
       solutionState=heatingSolution);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(StopTime=2000, __Dymola_Algorithm="Dassl"));
  end Definitions;

  model SimpleReaction
    extends Modelica.Icons.Example;
    Chemical.Processes.Reaction r(
      process=Chemical.Interfaces.processData(2),
      nS=1, nP=1) annotation (Placement(transformation(extent={{-6,14},{14,34}})));
    Chemical.Boundaries.Substance A(useFore=true) annotation (Placement(transformation(extent={{-50,14},{-30,34}})));
    Chemical.Boundaries.Substance B(useRear=true) annotation (Placement(transformation(extent={{30,14},{50,34}})));
  equation
    connect(A.fore, r.substrates[1]) annotation (Line(
        points={{-30,24},{-6,24}},
        color={158,66,200},
        thickness=0.5));
    connect(r.products[1], B.rear) annotation (Line(
        points={{14,24},{30,24}},
        color={158,66,200},
        thickness=0.5));
        annotation (Documentation(revisions="<html>
<p><i>2025</i></p>
<p>Marek Matejak</p>
</html>", info="<html>
        <p>Simple reaction demonstrating equilibration between substance A and substance B, in constant solution. Observe the molar concentration (A.c) and molar fraction.</p>
</html>"),
      experiment(StopTime=10));
  end SimpleReaction;

  model SimpleForwardReaction
    extends Modelica.Icons.Example;
    Processes.ForwardReaction forwardReaction(
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.Parameter,
      nS=1,
      nP=1) annotation (Placement(transformation(extent={{-8,14},{12,34}})));
    Chemical.Boundaries.Substance substance(
      useFore=true,
      use_mass_start=false,
      amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-70,14},{-50,34}})));
    Chemical.Boundaries.Substance substance1(
      useRear=true,
      use_mass_start=false,
      amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{48,12},{68,32}})));
    inner DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{-76,66},{-56,86}})));
  equation
    connect(substance.fore, forwardReaction.substrates[1]) annotation (Line(
        points={{-50,24},{-8,24}},
        color={158,66,200},
        thickness=0.5));
    connect(forwardReaction.products[1], substance1.rear)
      annotation (Line(
        points={{12,24},{40,24},{40,22},{48,22}},
        color={158,66,200},
        thickness=0.5));
    annotation (Documentation(revisions="<html>
<p><i>2025</i></p>
<p>Marek Matejak</p>
</html>", info="<html>
        <p>Simple forward reaction demonstrating one-way reaction from substance A to substance B, in constant solution. Observe the molar concentration (A.c) and molar fraction.</p>
</html>"),
      experiment(StopTime=10));
  end SimpleForwardReaction;

  model SimpleReactionInSolution "The simple chemical reaction A<->B with equilibrium B/A = 2"
    import Chemical;

     extends Modelica.Icons.Example;



    Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    Chemical.Boundaries.Substance A(
      useRear=false,
      useFore=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-52,-8},{-32,12}})));

    Chemical.Processes.Reaction
             reaction(
      process=Chemical.Interfaces.processData(2),

      nS=1,
      nP=1) annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
    Chemical.Boundaries.Substance B(
      useRear=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{44,-8},{64,12}})));

    inner Modelica.Fluid.System system annotation (Placement(transformation(extent={{58,64},{78,84}})));
  equation
    connect(A.solution, solution.solution) annotation (Line(
        points={{-48,-8},{-48,-86},{60,-86},{60,-98}},
        color={127,127,0}));
    connect(B.solution, solution.solution) annotation (Line(points={{48,-8},{48,-86},{60,-86},{60,-98}},
                                       color={127,127,0}));
    connect(A.fore, reaction.substrates[1]) annotation (Line(
        points={{-32,2},{-10,2}},
        color={158,66,200},
        thickness=0.5));
    connect(reaction.products[1], B.rear) annotation (Line(
        points={{10,2},{44,2}},
        color={158,66,200},
        thickness=0.5));
    annotation (Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Simple reaction demonstrating equilibria between substance A and substance B, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that mole fraction (A.x and B.x) are always summed to 1 for the solution.</p>
</html>"),
      experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
  end SimpleReactionInSolution;

  model SimpleReaction2 "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
    import Chemical;

     extends Modelica.Icons.Example;

    Chemical.Boundaries.Substance A(
      useFore=true,
      useSolution=false,
      use_mass_start=false,
      amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,4},{-14,24}})));
    Chemical.Processes.Reaction reaction2_1(
      process=Chemical.Interfaces.processData(2),
      nS=2,
      nP=1) annotation (Placement(transformation(extent={{6,-8},{26,12}})));
    Chemical.Boundaries.Substance B(
      useFore=true,
      useSolution=false,
      use_mass_start=false,
      amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
    Chemical.Boundaries.Substance C(
      useRear=true,
      useSolution=false,
      use_mass_start=false,
      amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{48,-8},{68,12}})));

  equation

    connect(A.fore, reaction2_1.substrates[1])
      annotation (Line(
        points={{-14,14},{-4,14},{-4,1.75},{6,1.75}},
        color={158,66,200},
        thickness=0.5));
    connect(B.fore, reaction2_1.substrates[2])
      annotation (Line(
        points={{-14,-14},{-8,-14},{-8,2.25},{6,2.25}},
        color={158,66,200},
        thickness=0.5));
    connect(reaction2_1.products[1], C.rear) annotation (Line(
        points={{26,2},{48,2}},
        color={158,66,200},
        thickness=0.5));
    annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
      experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
  end SimpleReaction2;

  model SimpleReaction22 "The simple chemical reaction A+B<->C+D with equilibrium [C]*[D]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
     extends Modelica.Icons.Example;

    Chemical.Boundaries.Substance A(
      useFore=true,                 use_mass_start=false, amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
    Chemical.Processes.Reaction reaction2_1(
      process=Chemical.Interfaces.processData(2),
      nS=2,
      nP=2)     annotation (Placement(transformation(extent={{2,-12},{22,8}})));
    Chemical.Boundaries.Substance B(
      useFore=true,                 use_mass_start=false, amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-32,-24},{-12,-4}})));
    Chemical.Boundaries.Substance C(useRear=true,
                                    amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{48,-8},{68,12}})));

    Chemical.Boundaries.Substance D(useRear=true,
                                    amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{44,-34},{64,-14}})));
    inner DropOfCommons dropOfCommons(L=1e-3) annotation (Placement(transformation(extent={{52,56},{72,76}})));
  equation

    connect(A.fore, reaction2_1.substrates[1])
      annotation (Line(
        points={{-14,12},{-4,12},{-4,-2.25},{2,-2.25}},
        color={158,66,200},
        thickness=0.5));
    connect(B.fore, reaction2_1.substrates[2])
      annotation (Line(
        points={{-12,-14},{-2,-14},{-2,-1.75},{2,-1.75}},
        color={158,66,200},
        thickness=0.5));
    connect(C.rear, reaction2_1.products[1]) annotation (Line(
        points={{48,2},{28,2},{28,-2.25},{22,-2.25}},
        color={158,66,200},
        thickness=0.5));
    connect(D.rear, reaction2_1.products[2])
      annotation (Line(
        points={{44,-24},{28,-24},{28,-1.75},{22,-1.75}},
        color={158,66,200},
        thickness=0.5));
    annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
      experiment(StopTime=100, __Dymola_Algorithm="Dassl"));
  end SimpleReaction22;

  model ExothermicReaction "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"
    import Chemical;

     extends Modelica.Icons.Example;

    parameter Modelica.Units.SI.MolarEnergy ReactionEnthalpy=-55000;

    Chemical.Solution thermal_isolated_solution(useMechanicPorts=true, ConstantTemperature=false)
      annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));
    Chemical.Boundaries.Substance A(
      useFore=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
    Chemical.Processes.Reaction reaction2_2(
      process=Chemical.Interfaces.processData(1, ReactionEnthalpy),
      nS=1,
      nP=1) annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
    Chemical.Boundaries.Substance B(
      useRear=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{28,-60},{48,-40}})));

    Chemical.Solution solution_at_constant_temperature(useMechanicPorts=true, useThermalPort=true)
      annotation (Placement(transformation(extent={{-100,0},{98,94}})));
    Chemical.Boundaries.Substance A1(
      useFore=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
    Chemical.Processes.Reaction reaction2_1(
      process=Chemical.Interfaces.processData(1, ReactionEnthalpy),
      nS=1,
      nP=1) annotation (Placement(transformation(extent={{-8,40},{12,60}})));
    Chemical.Boundaries.Substance B1(
      useRear=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{20,40},{40,60}})));

    //  Modelica.SIunits.HeatFlowRate q
    //    "Heat flow to environment to reach constant temperature";
    Modelica.Units.SI.Temperature t
      "Temperature if the solution is ideally thermal isolated from environment";
    Chemical.Boundaries.Substance H2O(substanceDefinition=Chemical.Substances.Liquid.H2O,
      useSolution=true,                                                                   mass_start=1)
      annotation (Placement(transformation(extent={{20,4},{40,24}})));
    Chemical.Boundaries.Substance H2O1(substanceDefinition=Chemical.Substances.Liquid.H2O,
      useSolution=true,                                                                    mass_start=1)
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

    connect(B.solution, thermal_isolated_solution.solution) annotation (Line(
        points={{32,-60},{32,-64},{58.4,-64},{58.4,-99.06}},
        color={127,127,0}));
    connect(A.solution, thermal_isolated_solution.solution) annotation (Line(
          points={{-36,-60},{-36,-64},{58.4,-64},{58.4,-99.06}},
                                                           color={127,127,0}));
    connect(B1.solution, solution_at_constant_temperature.solution) annotation (
        Line(
        points={{24,40},{24,34},{58.4,34},{58.4,0.94}},
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
    connect(A1.fore, reaction2_1.substrates[1]) annotation (Line(
        points={{-20,50},{-8,50}},
        color={158,66,200},
        thickness=0.5));
    connect(reaction2_1.products[1], B1.rear) annotation (Line(
        points={{12,50},{20,50}},
        color={158,66,200},
        thickness=0.5));
    connect(A.fore, reaction2_2.substrates[1]) annotation (Line(
        points={{-20,-50},{-8,-50}},
        color={158,66,200},
        thickness=0.5));
    connect(reaction2_2.products[1], B.rear) annotation (Line(
        points={{12,-50},{28,-50}},
        color={158,66,200},
        thickness=0.5));
    annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Demonstration of exotermic reaction with perfect cooling (i.e. connected fixed temperature to the HeatPort) and thermally insulated (HetPort unconnected). See solution_(...).T</p>
</html>"),
      experiment(StopTime=10, __Dymola_Algorithm="Dassl"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}})));
  end ExothermicReaction;

  model EnzymeKinetics "Basic enzyme kinetics"

    extends Modelica.Icons.Example;

    Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
    Chemical.Boundaries.Substance P(
      useRear=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=1e-8) annotation (Placement(transformation(extent={{72,-12},{92,8}})));

    Chemical.Boundaries.Substance S(
      useFore=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=100) annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));

    parameter Modelica.Units.SI.AmountOfSubstance tE=1
      "Total amount of enzyme";
       parameter Real k_cat(
      unit="mol/s",
      displayUnit="mol/min")=1
      "Forward rate of second reaction";
    constant Modelica.Units.SI.Concentration Km=0.1
      "Michaelis constant = substrate concentration at rate of half Vmax";

    parameter Modelica.Units.SI.MolarFlowRate Vmax=1e-5*k_cat
      "Maximal molar flow";

    Chemical.Boundaries.Substance ES(
      useRear=true,
      useFore=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=tE/10)
                                    annotation (Placement(transformation(extent={{-8,-10},{12,10}})));
    Chemical.Boundaries.Substance E(
      useRear=true,
      useFore=true,
      useSolution=true,
      use_mass_start=false,
      amountOfSubstance_start=tE*(9/10))
                                    annotation (Placement(transformation(extent={{12,36},{-8,56}})));
    Processes.Reaction chemicalReaction(
      process=Chemical.Interfaces.processData(2/Km),
      k_forward=1,
      nS=2,
      nP=1) annotation (Placement(transformation(extent={{-42,-8},{-22,12}})));

    Processes.ForwardReaction chemicalReaction1(
      k_forward=k_cat,
      nS=1,
      nP=2) annotation (Placement(transformation(extent={{26,-8},{46,12}})));

    Chemical.Boundaries.Substance liquidWater(
      useSolution=true,
      substanceDefinition=Chemical.Substances.Liquid.H2O,
      use_mass_start=true,
      mass_start=1) annotation (Placement(transformation(extent={{42,-80},{62,-60}})));
    inner DropOfCommons dropOfCommons      annotation (Placement(transformation(extent={{68,70},{88,90}})));
  equation
    //Michaelis-Menton: v=((E.q_out.conc + ES.q_out.conc)*k_cat)*S.concentration/(Km+S.concentration);
    connect(E.solution, solution.solution) annotation (Line(
        points={{8,36},{-8,36},{-8,-98},{60,-98}},
        color={127,127,0}));
    connect(ES.solution, solution.solution)
      annotation (Line(points={{-4,-10},{-4,-98},{60,-98}},         color={127,127,0}));

    connect(S.solution, solution.solution) annotation (Line(
        points={{-88,-14},{-88,-56},{-8,-56},{-8,-98},{60,-98}},
        color={127,127,0}));
    connect(P.solution, solution.solution) annotation (Line(
        points={{76,-12},{76,-98},{60,-98}},
        color={127,127,0}));
    connect(liquidWater.solution, solution.solution) annotation (Line(points={{
            46,-80},{46,-98},{60,-98}}, color={127,127,0}));
    connect(chemicalReaction.products[1], ES.rear) annotation (Line(
        points={{-22,2},{-16,2},{-16,0},{-8,0}},
        color={158,66,200},
        thickness=0.5));
    connect(ES.fore, chemicalReaction1.substrates[1])
      annotation (Line(
        points={{12,0},{18,0},{18,2},{26,2}},
        color={158,66,200},
        thickness=0.5));
    connect(S.fore, chemicalReaction.substrates[1])
      annotation (Line(
        points={{-72,-4},{-52,-4},{-52,-2},{-42,-2},{-42,1.75}},
        color={158,66,200},
        thickness=0.5));
    connect(E.fore, chemicalReaction.substrates[2])
      annotation (Line(
        points={{-8,46},{-52,46},{-52,2.25},{-42,2.25}},
        color={158,66,200},
        thickness=0.5));
    connect(chemicalReaction1.products[1], P.rear)
      annotation (Line(
        points={{46,1.75},{58,1.75},{58,-2},{72,-2}},
        color={158,66,200},
        thickness=0.5));
    connect(E.rear, chemicalReaction1.products[2])
      annotation (Line(
        points={{12,46},{48,46},{48,48},{56,48},{56,2.25},{46,2.25}},
        color={158,66,200},
        thickness=0.5));
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
      experiment(StopTime=100000, __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput);
  end EnzymeKinetics;

  model Diffusion

    extends Modelica.Icons.Example;

    Chemical.Boundaries.Substance s2(useFore=true, mass_start=0.6) annotation (Placement(transformation(extent={{-174,20},{-154,40}})));
    Chemical.Boundaries.Substance p2(useRear=true, mass_start=0.4) annotation (Placement(transformation(extent={{-66,20},{-46,40}})));
    Chemical.Processes.Diffusion d2(solutionFrom=Chemical.Utilities.Types.SolutionChoice.FirstSubstrate, redeclare function uLoss =
          Chemical.Processes.Internal.Kinetics.generalPotentialLoss) annotation (Placement(transformation(extent={{-118,20},{-98,40}})));
    Chemical.Boundaries.Substance s1(useFore=true, mass_start=0.6) annotation (Placement(transformation(extent={{-174,58},{-154,78}})));
    Chemical.Boundaries.Substance p1(useRear=true, mass_start=0.4) annotation (Placement(transformation(extent={{-66,58},{-46,78}})));
    Chemical.Processes.Diffusion d1(solutionFrom=Chemical.Utilities.Types.SolutionChoice.Parameter, redeclare function uLoss =
          Processes.Internal.Kinetics.generalPotentialLoss)
                                                  annotation (Placement(transformation(extent={{-118,58},{-98,78}})));
    Chemical.Boundaries.Substance s3(
      useFore=true,
      useSolution=true,
      mass_start=0.6) annotation (Placement(transformation(extent={{-170,-54},{-150,-34}})));
    Chemical.Boundaries.Substance p3(
      useRear=true,
      useSolution=true,
      mass_start=0.4) annotation (Placement(transformation(extent={{-62,-54},{-42,-34}})));
    Chemical.Processes.Diffusion d3(solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort, redeclare function uLoss =
          Processes.Internal.Kinetics.generalPotentialLoss)
                                                  annotation (Placement(transformation(extent={{-114,-56},{-94,-36}})));
    Solution solution annotation (Placement(transformation(extent={{-222,-122},{-14,-12}})));
    inner Modelica.Fluid.System system annotation (Placement(transformation(extent={{-210,64},{-190,84}})));
    Chemical.Boundaries.Substance s4(
      useFore=true,
      useSolution=false,
      mass_start=0.6) annotation (Placement(transformation(extent={{78,64},{98,84}})));
    Chemical.Boundaries.Substance p4(
      useRear=true,
      useSolution=false,
      mass_start=0.4) annotation (Placement(transformation(extent={{186,64},{206,84}})));
    Chemical.Processes.Diffusion d4(solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort, redeclare function uLoss =
          Processes.Internal.Kinetics.generalPotentialLoss)
                                                  annotation (Placement(transformation(extent={{134,62},{154,82}})));
    Solution solution1 annotation (Placement(transformation(extent={{26,-4},{234,106}})));
    Chemical.Boundaries.Substance solvent(useFore=false, useSolution=true,
      mass_start=1)                                                        annotation (Placement(transformation(extent={{194,24},{214,44}})));
    Chemical.Boundaries.Substance s5(
      useFore=true,
      useSolution=true,
      mass_start=0.6) annotation (Placement(transformation(extent={{64,-54},{84,-34}})));
    Chemical.Boundaries.Substance p5(
      useRear=true,
      useSolution=true,
      mass_start=0.4) annotation (Placement(transformation(extent={{172,-54},{192,-34}})));
    Chemical.Processes.Diffusion d5(solutionFrom=Chemical.Utilities.Types.SolutionChoice.FirstSubstrate, redeclare function uLoss =
          Processes.Internal.Kinetics.generalPotentialLoss)
                                                  annotation (Placement(transformation(extent={{120,-56},{140,-36}})));
    Solution solution2 annotation (Placement(transformation(extent={{12,-122},{220,-12}})));
  equation
    connect(s2.fore, d2.rear) annotation (Line(
        points={{-154,30},{-118,30}},
        color={158,66,200},
        thickness=0.5));
    connect(d2.fore, p2.rear) annotation (Line(
        points={{-98,30},{-66,30}},
        color={158,66,200},
        thickness=0.5));
    connect(s1.fore, d1.rear) annotation (Line(
        points={{-154,68},{-118,68}},
        color={158,66,200},
        thickness=0.5));
    connect(d1.fore, p1.rear) annotation (Line(
        points={{-98,68},{-66,68}},
        color={158,66,200},
        thickness=0.5));
    connect(s3.fore, d3.rear) annotation (Line(
        points={{-150,-44},{-122,-44},{-122,-46},{-114,-46}},
        color={158,66,200},
        thickness=0.5));
    connect(d3.fore, p3.rear) annotation (Line(
        points={{-94,-46},{-70,-46},{-70,-44},{-62,-44}},
        color={158,66,200},
        thickness=0.5));
    connect(solution.solution, d3.solution) annotation (Line(points={{-55.6,-120.9},{-55.6,-126},{-110,-126},{-110,-56}}, color={127,127,0}));
    connect(s3.solution, solution.solution) annotation (Line(points={{-166,-54},{-166,-120.9},{-55.6,-120.9}}, color={127,127,0}));
    connect(p3.solution, solution.solution) annotation (Line(points={{-58,-54},{-58,-126},{-55.6,-126},{-55.6,-120.9}}, color={127,127,0}));
    connect(s4.fore, d4.rear) annotation (Line(
        points={{98,74},{126,74},{126,72},{134,72}},
        color={158,66,200},
        thickness=0.5));
    connect(d4.fore, p4.rear) annotation (Line(
        points={{154,72},{178,72},{178,74},{186,74}},
        color={158,66,200},
        thickness=0.5));
    connect(solution1.solution, d4.solution) annotation (Line(points={{192.4,-2.9},{192.4,-8},{138,-8},{138,62}}, color={127,127,0}));
    connect(solution1.solution, solvent.solution) annotation (Line(points={{192.4,-2.9},{192.4,-10},{198,-10},{198,24}},     color={127,127,0}));
    connect(s5.fore, d5.rear) annotation (Line(
        points={{84,-44},{112,-44},{112,-46},{120,-46}},
        color={158,66,200},
        thickness=0.5));
    connect(d5.fore, p5.rear) annotation (Line(
        points={{140,-46},{164,-46},{164,-44},{172,-44}},
        color={158,66,200},
        thickness=0.5));
    connect(s5.solution, solution2.solution) annotation (Line(points={{68,-54},{68,-120.9},{178.4,-120.9}}, color={127,127,0}));
    connect(p5.solution, solution2.solution) annotation (Line(points={{176,-54},{176,-126},{178.4,-126},{178.4,-120.9}}, color={127,127,0}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-240,-140},{240,140}})), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-240,
              -140},{240,140}})),
      experiment(StopTime=1, __Dymola_Algorithm="Dassl"));
  end Diffusion;

  model SimpleFlow
    extends Modelica.Icons.Example;
    constant Real K = 2 "Dissociation constant of the reaction";

    constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
    constant Real R = Modelica.Constants.R "Gas constant";
    Chemical.Boundaries.Substance substance(
      useFore=true,
      use_mass_start=false,
      amountOfSubstance_start=2)   annotation (Placement(transformation(extent={{-70,14},{-50,34}})));
    Chemical.Boundaries.Substance substance1(
      useRear=true,
      use_mass_start=false,
      amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{48,12},{68,32}})));
    inner DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{-76,66},{-56,86}})));
    Processes.Pump pump(solutionFrom=Chemical.Utilities.Types.SolutionChoice.FirstSubstrate, SubstanceFlow=1)
      annotation (Placement(transformation(extent={{-14,14},{6,34}})));
  equation
    connect(substance.fore, pump.rear) annotation (Line(
        points={{-50,24},{-14,24}},
        color={158,66,200},
        thickness=0.5));
    connect(pump.fore, substance1.rear) annotation (Line(
        points={{6,24},{40,24},{40,22},{48,22}},
        color={158,66,200},
        thickness=0.5));
  end SimpleFlow;

  model SimpleReactionPathway
    extends Modelica.Icons.Example;

    Chemical.Processes.Reaction r1(
      process=Chemical.Interfaces.processData(2),
      nS=1,
      nP=1) annotation (Placement(transformation(extent={{-42,48},{-22,68}})));
    Chemical.Boundaries.Substance A(useFore=true) annotation (Placement(transformation(extent={{-76,48},{-56,68}})));
    Chemical.Boundaries.Substance B(useRear=true) annotation (Placement(transformation(extent={{60,48},{80,68}})));
    Processes.Reaction r2(nP=1, nS=1) annotation (Placement(transformation(extent={{-4,48},{16,68}})));
    Processes.GasSolubility gasSolubility annotation (Placement(transformation(extent={{28,48},{48,68}})));
  equation
    connect(A.fore, r1.substrates[1]) annotation (Line(
        points={{-56,58},{-42,58}},
        color={158,66,200},
        thickness=0.5));
    connect(r2.products[1], gasSolubility.rear) annotation (Line(
        points={{16,58},{28,58}},
        color={158,66,200},
        thickness=0.5));
    connect(gasSolubility.fore, B.rear) annotation (Line(
        points={{48,58},{60,58}},
        color={158,66,200},
        thickness=0.5));
    connect(r1.products[1], r2.substrates[1]) annotation (Line(
        points={{-22,58},{-4,58}},
        color={158,66,200},
        thickness=0.5));
  end SimpleReactionPathway;
  annotation (Documentation(info="<html>
<u>Tests for top level components of the undirected chemical simulation package.</u>
</html>"));
end Examples;
