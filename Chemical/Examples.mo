within Chemical;
package Examples "Tests for top level components of undirected"
  extends Modelica.Icons.ExamplesPackage;

  model Definitions
    extends Modelica.Icons.Example;
    import Chemical.Interfaces.Properties;
    import Chemical.Interfaces.Definition;
    import Chemical.Substances.Gas;
    import Chemical.Substances.Liquid;
    import Chemical.Interfaces.processData;
    constant Real R = Modelica.Constants.R;
    constant Definition O2_aq = Gas.O2 + processData(0.0013, -1500*R);
    constant Definition NO_aq = Gas.NO + processData(0.0014, -1500*R);
    constant Definition H2O_formation = Gas.H2O - (Gas.H2 + 0.5*Gas.O2);
    constant Definition Hemoglobin = Liquid.Unknown;
    import Chemical.Interfaces.SolutionState;
    import Chemical.Interfaces.Phase;
    constant SolutionState SATP = SolutionState(phase = Phase.Gas, T = 298.15);
    SolutionState heatingSolution = SolutionState(phase = Phase.Gas, T = 273.15+1*time);
    import Chemical.Interfaces.ProcessProperties;
    ProcessProperties O2_dissolving_props(definition = processData(0.0013, -1500*R), solutionState = heatingSolution);
    ProcessProperties H2O_formation_props(definition = H2O_formation, solutionState = heatingSolution);
    ProcessProperties H2O_vaporization_props(definition = Gas.H2O - Liquid.H2O, solutionState = heatingSolution);
    import Modelica.Units.SI;
    SI.ChemicalPotential uO2,uH2,uH2O,uNO,uO2_,uNO_,uNO_aq;
    SI.MolarEnthalpy hO2;
    SI.MolarEntropy sO2;
  equation
    uO2 = Properties.electroChemicalPotentialPure(Gas.O2,heatingSolution);
    uNO = Properties.electroChemicalPotentialPure(Gas.NO,heatingSolution);
    uO2_ = Properties.electroChemicalPotentialPure(Gas.O2,SATP);
    uNO_ = Properties.electroChemicalPotentialPure(Gas.NO,SATP);
    uNO_aq = Properties.electroChemicalPotentialPure(NO_aq,SATP);

    hO2 = Properties.molarEnthalpyElectroneutral(Gas.O2,heatingSolution);
    sO2 = Properties.molarEntropyPure(Gas.O2,heatingSolution);
    uH2 = Chemical.Interfaces.Properties.electroChemicalPotentialPure(Gas.H2,heatingSolution);
    uH2O = Chemical.Interfaces.Properties.electroChemicalPotentialPure(Gas.H2O,heatingSolution);
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false)),
      Diagram(coordinateSystem(preserveAspectRatio = false)),
      experiment(StopTime = 2000, __Dymola_Algorithm = "Dassl"));
  end Definitions;

  model Diffusion
    extends Modelica.Icons.Example;
    Chemical.Boundaries.Substance s2(useFore = true, mass_start = 0.6) annotation(
      Placement(transformation(extent = {{-174, 20}, {-154, 40}})));
    Chemical.Boundaries.Substance p2(useRear = true, mass_start = 0.4) annotation(
      Placement(transformation(extent = {{-66, 20}, {-46, 40}})));
    Chemical.Processes.Diffusion d2(solutionFrom = Chemical.Utilities.Types.SolutionChoice.FirstSubstrate, redeclare function uLoss = Chemical.Processes.Internal.Kinetics.generalPotentialLoss) annotation(
      Placement(transformation(extent = {{-118, 20}, {-98, 40}})));
    Chemical.Boundaries.Substance s1(useFore = true, mass_start = 0.6) annotation(
      Placement(transformation(extent = {{-174, 58}, {-154, 78}})));
    Chemical.Boundaries.Substance p1(useRear = true, mass_start = 0.4) annotation(
      Placement(transformation(extent = {{-66, 58}, {-46, 78}})));
    Chemical.Processes.Diffusion d1(solutionFrom = Chemical.Utilities.Types.SolutionChoice.Parameter, redeclare function uLoss = Processes.Internal.Kinetics.generalPotentialLoss) annotation(
      Placement(transformation(extent = {{-118, 58}, {-98, 78}})));
    Chemical.Boundaries.Substance s3(useFore = true, useSolution = true, mass_start = 0.6) annotation(
      Placement(transformation(extent = {{-170, -54}, {-150, -34}})));
    Chemical.Boundaries.Substance p3(useRear = true, useSolution = true, mass_start = 0.4) annotation(
      Placement(transformation(extent = {{-62, -54}, {-42, -34}})));
    Chemical.Processes.Diffusion d3(solutionFrom = Chemical.Utilities.Types.SolutionChoice.SolutionPort, redeclare function uLoss = Processes.Internal.Kinetics.generalPotentialLoss) annotation(
      Placement(transformation(extent = {{-114, -56}, {-94, -36}})));
    Solution solution annotation(
      Placement(transformation(extent = {{-222, -122}, {-14, -12}})));
    inner Modelica.Fluid.System system annotation(
      Placement(transformation(extent = {{-210, 64}, {-190, 84}})));
    Chemical.Boundaries.Substance s4(useFore = true, useSolution = false, mass_start = 0.6) annotation(
      Placement(transformation(extent = {{78, 64}, {98, 84}})));
    Chemical.Boundaries.Substance p4(useRear = true, useSolution = false, mass_start = 0.4) annotation(
      Placement(transformation(extent = {{186, 64}, {206, 84}})));
    Chemical.Processes.Diffusion d4(solutionFrom = Chemical.Utilities.Types.SolutionChoice.SolutionPort, redeclare function uLoss = Processes.Internal.Kinetics.generalPotentialLoss) annotation(
      Placement(transformation(extent = {{134, 62}, {154, 82}})));
    Solution solution1 annotation(
      Placement(transformation(extent = {{26, -4}, {234, 106}})));
    Chemical.Boundaries.Substance solvent(useFore = false, useSolution = true, mass_start = 1) annotation(
      Placement(transformation(extent = {{194, 24}, {214, 44}})));
    Chemical.Boundaries.Substance s5(useFore = true, useSolution = true, mass_start = 0.6) annotation(
      Placement(transformation(extent = {{64, -54}, {84, -34}})));
    Chemical.Boundaries.Substance p5(useRear = true, useSolution = true, mass_start = 0.4) annotation(
      Placement(transformation(extent = {{172, -54}, {192, -34}})));
    Chemical.Processes.Diffusion d5(solutionFrom = Chemical.Utilities.Types.SolutionChoice.FirstSubstrate, redeclare function uLoss = Processes.Internal.Kinetics.generalPotentialLoss) annotation(
      Placement(transformation(extent = {{120, -56}, {140, -36}})));
    Solution solution2 annotation(
      Placement(transformation(extent = {{12, -122}, {220, -12}})));
  equation
    connect(s2.fore, d2.rear) annotation(
      Line(points = {{-154, 30}, {-118, 30}}, color = {158, 66, 200}, thickness = 0.5));
    connect(d2.fore, p2.rear) annotation(
      Line(points = {{-98, 30}, {-66, 30}}, color = {158, 66, 200}, thickness = 0.5));
    connect(s1.fore, d1.rear) annotation(
      Line(points = {{-154, 68}, {-118, 68}}, color = {158, 66, 200}, thickness = 0.5));
    connect(d1.fore, p1.rear) annotation(
      Line(points = {{-98, 68}, {-66, 68}}, color = {158, 66, 200}, thickness = 0.5));
    connect(s3.fore, d3.rear) annotation(
      Line(points = {{-150, -44}, {-122, -44}, {-122, -46}, {-114, -46}}, color = {158, 66, 200}, thickness = 0.5));
    connect(d3.fore, p3.rear) annotation(
      Line(points = {{-94, -46}, {-70, -46}, {-70, -44}, {-62, -44}}, color = {158, 66, 200}, thickness = 0.5));
    connect(solution.solution, d3.solution) annotation(
      Line(points = {{-55.6, -120.9}, {-55.6, -126}, {-110, -126}, {-110, -56}}, color = {127, 127, 0}));
    connect(s3.solution, solution.solution) annotation(
      Line(points = {{-166, -54}, {-166, -120.9}, {-55.6, -120.9}}, color = {127, 127, 0}));
    connect(p3.solution, solution.solution) annotation(
      Line(points = {{-58, -54}, {-58, -126}, {-55.6, -126}, {-55.6, -120.9}}, color = {127, 127, 0}));
    connect(s4.fore, d4.rear) annotation(
      Line(points = {{98, 74}, {126, 74}, {126, 72}, {134, 72}}, color = {158, 66, 200}, thickness = 0.5));
    connect(d4.fore, p4.rear) annotation(
      Line(points = {{154, 72}, {178, 72}, {178, 74}, {186, 74}}, color = {158, 66, 200}, thickness = 0.5));
    connect(solution1.solution, d4.solution) annotation(
      Line(points = {{192.4, -2.9}, {192.4, -8}, {138, -8}, {138, 62}}, color = {127, 127, 0}));
    connect(solution1.solution, solvent.solution) annotation(
      Line(points = {{192.4, -2.9}, {192.4, -10}, {198, -10}, {198, 24}}, color = {127, 127, 0}));
    connect(s5.fore, d5.rear) annotation(
      Line(points = {{84, -44}, {112, -44}, {112, -46}, {120, -46}}, color = {158, 66, 200}, thickness = 0.5));
    connect(d5.fore, p5.rear) annotation(
      Line(points = {{140, -46}, {164, -46}, {164, -44}, {172, -44}}, color = {158, 66, 200}, thickness = 0.5));
    connect(s5.solution, solution2.solution) annotation(
      Line(points = {{68, -54}, {68, -120.9}, {178.4, -120.9}}, color = {127, 127, 0}));
    connect(p5.solution, solution2.solution) annotation(
      Line(points = {{176, -54}, {176, -126}, {178.4, -126}, {178.4, -120.9}}, color = {127, 127, 0}));
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-240, -140}, {240, 140}})),
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-240, -140}, {240, 140}})),
      experiment(StopTime = 1, __Dymola_Algorithm = "Dassl"));
  end Diffusion;

  model OxygenInWater "Oxygen dsisolution in water"
    extends Modelica.Icons.Example;


    Boundaries.Substance O2aqS(
      useRear=true,
      useFore=false,
      useSolution=false,
      preferMass=false,
      amountOfSubstance_start=2e-4) annotation (Placement(transformation(extent={{60,2},{80,22}})));
    Boundaries.Substance O2gS(
      substanceDefinition=Chemical.Substances.Gas.O2,
      useFore=true,
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Gas, T=310.15),
      preferMass=false,
      amountOfSubstance_start=0.2*1e5*1e-3/(Modelica.Constants.R*(273.15 + 37))) annotation (Placement(transformation(extent={{-38,2},{-18,22}})));
    Boundaries.Substance O2aqR(
      useRear=true,
      useFore=false,
      useSolution=false,
      preferMass=false,
      amountOfSubstance_start=2e-4) annotation (Placement(transformation(extent={{68,-30},{88,-10}})));
    Boundaries.Substance O2gR(
      substanceDefinition=Chemical.Substances.Gas.O2,
      useFore=true,
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Gas, T=310.15),
      preferMass=false,
      amountOfSubstance_start=0.2*1e5*1e-3/(Modelica.Constants.R*(273.15 + 37))) annotation (Placement(transformation(extent={{-36,-30},{-16,-10}})));
    Processes.Reaction gasSolubilityR(
      firstProductFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
      firstProduct=Chemical.Substances.Aqueous.O2,
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.Parameter,
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible, T=310.15),
      nS=1,
      nP=1) annotation (Placement(transformation(extent={{16,-30},{36,-10}})));
    Processes.GasSolubility gasSolubilityS(
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible, T=310.15),
      productFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
      product=Chemical.Substances.Aqueous.O2) annotation (Placement(transformation(extent={{12,2},{32,22}})));
    Boundaries.ExternalGas O2gE(substanceDefinition=Chemical.Substances.Gas.O2, PartialPressure(displayUnit="mmHg") = 13332.2387415)
      annotation (Placement(transformation(extent={{-36,66},{-16,86}})));
    Boundaries.Substance O2aqE(
      useRear=true,
      useFore=false,
      useSolution=false,
      preferMass=false,
      amountOfSubstance_start=2e-4) annotation (Placement(transformation(extent={{56,66},{76,86}})));
    Processes.GasSolubility gasSolubilityE(
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible, T=310.15),
      productFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
      product=Chemical.Substances.Aqueous.O2) annotation (Placement(transformation(extent={{12,66},{32,86}})));
    Boundaries.Substance O2aqF(
      useRear=true,
      useFore=false,
      useSolution=false,
      preferMass=false,
      amountOfSubstance_start=2e-9) annotation (Placement(transformation(extent={{70,-96},{90,-76}})));
    Processes.GasSolubility gasSolubilityF(
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible, T=310.15),
      productFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
      product=Chemical.Substances.Aqueous.O2) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={28,-86})));
    Boundaries.TerminalInflow O2gF(
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Gas, T=310.15),
      SubstanceFlow(displayUnit="umol/min") = 1.6666666666667e-10,
      definition=Chemical.Substances.Gas.O2) annotation (Placement(transformation(extent={{-26,-96},{-6,-76}})));
    Modelica.Blocks.Sources.ContinuousClock clock(offset=1e3)
                                                            annotation (
       Placement(transformation(extent={{-58,-62},{-38,-42}})));
    Boundaries.Substance O2aqT(
      useRear=true,
      useFore=false,
      useSolution=false,
      preferMass=false,
      amountOfSubstance_start=2e-3) annotation (Placement(transformation(extent={{58,-70},{80,-50}})));
    Boundaries.ExternalGas O2gT(
      substanceDefinition=Chemical.Substances.Gas.O2,
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Gas, T=310.15),
      usePartialPressureInput=true,
      PartialPressure(displayUnit="mmHg") = 13332.2387415)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-14,-60})));
    Processes.GasSolubility gasSolubilityT(
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible, T=310.15),
      productFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
      product=Chemical.Substances.Aqueous.O2) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={28,-60})));
    Boundaries.ExternalGas O2gP(substanceDefinition=Chemical.Substances.Gas.O2, PartialPressure(displayUnit="mmHg") = 13332.2387415)
      annotation (Placement(transformation(extent={{-36,32},{-16,52}})));
    Boundaries.Substance O2aqP(
      useRear=true,
      useFore=false,
      useSolution=false,
      preferMass=false,
      amountOfSubstance_start=2e-4) annotation (Placement(transformation(extent={{56,32},{76,52}})));
    Processes.GasSolubility gasSolubilityP(
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible, T=310.15),
      productFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
      product=Chemical.Substances.Aqueous.O2,
      process=Chemical.Interfaces.processData(0.0013, -1500*Modelica.Constants.R))
                                              annotation (Placement(transformation(extent={{12,32},{32,52}})));
  equation

    connect(O2gR.fore, gasSolubilityR.substrates[1]) annotation (Line(
        points={{-16,-20},{16,-20}},
        color={158,66,200},
        thickness=0.5));
    connect(gasSolubilityR.products[1], O2aqR.rear) annotation (Line(
        points={{36,-20},{68,-20}},
        color={158,66,200},
        thickness=0.5));
    connect(O2gS.fore, gasSolubilityS.rear) annotation (Line(
        points={{-18,12},{12,12}},
        color={158,66,200},
        thickness=0.5));
    connect(gasSolubilityS.fore, O2aqS.rear) annotation (Line(
        points={{32,12},{60,12}},
        color={158,66,200},
        thickness=0.5));
    connect(gasSolubilityE.fore, O2aqE.rear) annotation (Line(
        points={{32,76},{56,76}},
        color={158,66,200},
        thickness=0.5));
    connect(O2gE.fore,gasSolubilityE. rear) annotation (Line(
        points={{-16,76},{12,76}},
        color={158,66,200},
        thickness=0.5));
    connect(gasSolubilityF.fore, O2aqF.rear) annotation (Line(
        points={{38,-86},{70,-86}},
        color={158,66,200},
        thickness=0.5));
    connect(O2gF.fore, gasSolubilityF.rear) annotation (Line(
        points={{-6,-86},{18,-86}},
        color={158,66,200},
        thickness=0.5));
    connect(gasSolubilityT.rear, O2gT.fore) annotation (Line(
        points={{18,-60},{-4,-60}},
        color={158,66,200},
        thickness=0.5));
    connect(gasSolubilityT.fore, O2aqT.rear) annotation (Line(
        points={{38,-60},{58,-60}},
        color={158,66,200},
        thickness=0.5));
    connect(clock.y, O2gT.partialPressure) annotation (Line(points={{-37,-52},{-36,-52.8},{-24,-52.8}}, color={0,0,127}));
    connect(gasSolubilityP.fore, O2aqP.rear) annotation (Line(
        points={{32,42},{56,42}},
        color={158,66,200},
        thickness=0.5));
    connect(O2gP.fore, gasSolubilityP.rear) annotation (Line(
        points={{-16,42},{12,42}},
        color={158,66,200},
        thickness=0.5));
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
  end OxygenInWater;

  model SimpleFlow
    extends Modelica.Icons.Example;
    constant Real K = 2 "Dissociation constant of the reaction";
    constant Modelica.Units.SI.Temperature T_25degC = 298.15 "Temperature";
    constant Real R = Modelica.Constants.R "Gas constant";
    Chemical.Boundaries.Substance substance(useFore = true, preferMass = false, amountOfSubstance_start = 2) annotation(
      Placement(transformation(extent = {{-70, 14}, {-50, 34}})));
    Chemical.Boundaries.Substance substance1(useRear = true, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{48, 12}, {68, 32}})));
    inner DropOfCommons dropOfCommons annotation(
      Placement(transformation(extent = {{-76, 66}, {-56, 86}})));
    Processes.Pump pump(solutionFrom = Chemical.Utilities.Types.SolutionChoice.FirstSubstrate, SubstanceFlow = 1) annotation(
      Placement(transformation(extent = {{-14, 14}, {6, 34}})));
  equation
    connect(substance.fore, pump.rear) annotation(
      Line(points = {{-50, 24}, {-14, 24}}, color = {158, 66, 200}, thickness = 0.5));
    connect(pump.fore, substance1.rear) annotation(
      Line(points = {{6, 24}, {40, 24}, {40, 22}, {48, 22}}, color = {158, 66, 200}, thickness = 0.5));
  end SimpleFlow;

  model SimpleReaction
    extends Modelica.Icons.Example;
    Chemical.Processes.Reaction r( nP = 1, nS = 1, process = Chemical.Interfaces.processData(2)) annotation(
      Placement(transformation(extent={{-6,12},{14,32}})));
    Chemical.Boundaries.Substance A(useFore = true) annotation(
      Placement(transformation(extent = {{-50, 14}, {-30, 34}})));
    Chemical.Boundaries.Substance B(useRear = true) annotation(
      Placement(transformation(extent = {{30, 14}, {50, 34}})));
  equation
    connect(A.fore, r.substrates[1]) annotation(
      Line(points={{-30,24},{-18,24},{-18,22},{-6,22}},
                                           color = {158, 66, 200}, thickness = 0.5));
    connect(r.products[1], B.rear) annotation(
      Line(points={{14,22},{22,22},{22,24},{30,24}},
                                          color = {158, 66, 200}, thickness = 0.5));
    annotation(
      Documentation(revisions = "<html>
<p><i>2025</i></p>
<p>Marek Matejak</p>
</html>", info = "<html>
        <p>Simple reaction demonstrating equilibration between substance A and substance B, in constant solution. Observe the molar concentration (A.c) and molar fraction.</p>
</html>"),
      experiment(StopTime = 10));
  end SimpleReaction;

  model SimpleForwardReaction
    extends Modelica.Icons.Example;
    Processes.ForwardReaction forwardReaction(solutionFrom = Chemical.Utilities.Types.SolutionChoice.Parameter, nS = 1, nP = 1) annotation(
      Placement(transformation(extent = {{-8, 14}, {12, 34}})));
    Chemical.Boundaries.Substance substance(useFore = true, preferMass = false, amountOfSubstance_start = 0.9) annotation(
      Placement(transformation(extent = {{-70, 14}, {-50, 34}})));
    Chemical.Boundaries.Substance substance1(useRear = true, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{48, 12}, {68, 32}})));
    inner DropOfCommons dropOfCommons annotation(
      Placement(transformation(extent = {{-76, 66}, {-56, 86}})));
  equation
    connect(substance.fore, forwardReaction.substrates[1]) annotation(
      Line(points = {{-50, 24}, {-8, 24}}, color = {158, 66, 200}, thickness = 0.5));
    connect(forwardReaction.products[1], substance1.rear) annotation(
      Line(points = {{12, 24}, {40, 24}, {40, 22}, {48, 22}}, color = {158, 66, 200}, thickness = 0.5));
    annotation(
      Documentation(revisions = "<html>
<p><i>2025</i></p>
<p>Marek Matejak</p>
</html>", info = "<html>
        <p>Simple forward reaction demonstrating one-way reaction from substance A to substance B, in constant solution. Observe the molar concentration (A.c) and molar fraction.</p>
</html>"),
      experiment(StopTime = 10));
  end SimpleForwardReaction;

  model SimpleReactionPathway
    extends Modelica.Icons.Example;
    Chemical.Processes.Reaction r1(process = Chemical.Interfaces.processData(2), nS = 1, nP = 1) annotation(
      Placement(transformation(extent = {{-42, 48}, {-22, 68}})));
    Chemical.Boundaries.Substance A(useFore = true) annotation(
      Placement(transformation(extent = {{-76, 48}, {-56, 68}})));
    Chemical.Boundaries.Substance B(useRear = true) annotation(
      Placement(transformation(extent = {{60, 48}, {80, 68}})));
    Processes.Reaction r2(nP = 1, nS = 1) annotation(
      Placement(transformation(extent = {{-4, 48}, {16, 68}})));
    Processes.GasSolubility gasSolubility annotation(
      Placement(transformation(extent = {{28, 48}, {48, 68}})));
  equation
    connect(A.fore, r1.substrates[1]) annotation(
      Line(points = {{-56, 58}, {-42, 58}}, color = {158, 66, 200}, thickness = 0.5));
    connect(r2.products[1], gasSolubility.rear) annotation(
      Line(points = {{16, 58}, {28, 58}}, color = {158, 66, 200}, thickness = 0.5));
    connect(gasSolubility.fore, B.rear) annotation(
      Line(points = {{48, 58}, {60, 58}}, color = {158, 66, 200}, thickness = 0.5));
    connect(r1.products[1], r2.substrates[1]) annotation(
      Line(points = {{-22, 58}, {-4, 58}}, color = {158, 66, 200}, thickness = 0.5));
  end SimpleReactionPathway;

  model SimpleReactionInSolution "The simple chemical reaction A<->B with equilibrium B/A = 2"
    import Chemical;
    extends Modelica.Icons.Example;
    Chemical.Solution solution annotation(
      Placement(transformation(extent = {{-100, -100}, {100, 100}})));
    Chemical.Boundaries.Substance A(useRear = false, useFore = true, useSolution = true, preferMass = false, amountOfSubstance_start = 0.9) annotation(
      Placement(transformation(extent = {{-52, -8}, {-32, 12}})));
    Chemical.Processes.Reaction reaction(process = Chemical.Interfaces.processData(2), nS = 1, nP = 1) annotation(
      Placement(transformation(extent = {{-10, -8}, {10, 12}})));
    Chemical.Boundaries.Substance B(useRear = true, useSolution = true, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{44, -8}, {64, 12}})));
    inner Modelica.Fluid.System system annotation(
      Placement(transformation(extent = {{58, 64}, {78, 84}})));
  equation
    connect(A.solution, solution.solution) annotation(
      Line(points = {{-48, -8}, {-48, -86}, {60, -86}, {60, -98}}, color = {127, 127, 0}));
    connect(B.solution, solution.solution) annotation(
      Line(points = {{48, -8}, {48, -86}, {60, -86}, {60, -98}}, color = {127, 127, 0}));
    connect(A.fore, reaction.substrates[1]) annotation(
      Line(points = {{-32, 2}, {-10, 2}}, color = {158, 66, 200}, thickness = 0.5));
    connect(reaction.products[1], B.rear) annotation(
      Line(points = {{10, 2}, {44, 2}}, color = {158, 66, 200}, thickness = 0.5));
    annotation(
      Documentation(revisions = "<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info = "<html>
<p>Simple reaction demonstrating equilibria between substance A and substance B, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that mole fraction (A.x and B.x) are always summed to 1 for the solution.</p>
</html>"),
      experiment(StopTime = 10, __Dymola_Algorithm = "Dassl"));
  end SimpleReactionInSolution;

  model SimpleReaction2 "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
    import Chemical;
    extends Modelica.Icons.Example;
    Chemical.Boundaries.Substance A(useFore = true, useSolution = false, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{-34, 4}, {-14, 24}})));
    Chemical.Processes.Reaction reaction2_1(process = Chemical.Interfaces.processData(2), nS = 2, nP = 1) annotation(
      Placement(transformation(extent = {{6, -8}, {26, 12}})));
    Chemical.Boundaries.Substance B(useFore = true, useSolution = false, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{-34, -24}, {-14, -4}})));
    Chemical.Boundaries.Substance C(useRear = true, useSolution = false, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{48, -8}, {68, 12}})));
  equation
    connect(A.fore, reaction2_1.substrates[1]) annotation(
      Line(points = {{-14, 14}, {-4, 14}, {-4, 1.75}, {6, 1.75}}, color = {158, 66, 200}, thickness = 0.5));
    connect(B.fore, reaction2_1.substrates[2]) annotation(
      Line(points = {{-14, -14}, {-8, -14}, {-8, 2.25}, {6, 2.25}}, color = {158, 66, 200}, thickness = 0.5));
    connect(reaction2_1.products[1], C.rear) annotation(
      Line(points = {{26, 2}, {48, 2}}, color = {158, 66, 200}, thickness = 0.5));
    annotation(
      Documentation(revisions = "<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info = "<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
      experiment(StopTime = 10, __Dymola_Algorithm = "Dassl"));
  end SimpleReaction2;

  model SimpleReaction22 "The simple chemical reaction A+B<->C+D with equilibrium [C]*[D]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
    extends Modelica.Icons.Example;
    Chemical.Boundaries.Substance A(useFore = true, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{-34, 2}, {-14, 22}})));
    Chemical.Processes.Reaction reaction2_1(process = Chemical.Interfaces.processData(2), nS = 2, nP = 2) annotation(
      Placement(transformation(extent = {{2, -12}, {22, 8}})));
    Chemical.Boundaries.Substance B(useFore = true, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{-32, -24}, {-12, -4}})));
    Chemical.Boundaries.Substance C(useRear = true, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{48, -8}, {68, 12}})));
    Chemical.Boundaries.Substance D(useRear = true, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{44, -34}, {64, -14}})));
    inner DropOfCommons dropOfCommons           annotation(
      Placement(transformation(extent={{56,56},{76,76}})));
  equation
    connect(A.fore, reaction2_1.substrates[1]) annotation(
      Line(points = {{-14, 12}, {-4, 12}, {-4, -2.25}, {2, -2.25}}, color = {158, 66, 200}, thickness = 0.5));
    connect(B.fore, reaction2_1.substrates[2]) annotation(
      Line(points = {{-12, -14}, {-2, -14}, {-2, -1.75}, {2, -1.75}}, color = {158, 66, 200}, thickness = 0.5));
    connect(C.rear, reaction2_1.products[1]) annotation(
      Line(points = {{48, 2}, {28, 2}, {28, -2.25}, {22, -2.25}}, color = {158, 66, 200}, thickness = 0.5));
    connect(D.rear, reaction2_1.products[2]) annotation(
      Line(points = {{44, -24}, {28, -24}, {28, -1.75}, {22, -1.75}}, color = {158, 66, 200}, thickness = 0.5));
    annotation(
      Documentation(revisions = "<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info = "<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
      experiment(StopTime=0.0002, __Dymola_Algorithm="Dassl"));
  end SimpleReaction22;

  model HeatingOfWater "Heating of 1 kg water"
    import Chemical;
    extends Modelica.Icons.Example;

    Chemical.Solution solution(useMechanicPorts=true, useThermalPort=true) annotation (Placement(transformation(extent={{-98,-100},{102,100}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
      annotation (Placement(transformation(extent={{-86,-72},{-66,-52}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
    Chemical.Boundaries.Substance liquidWater(
      substanceDefinition=Chemical.Substances.Liquid.H2O,
      useFore=false,
      useSolution=true,
      preferMass=true,
      mass_start=1) annotation (Placement(transformation(extent={{22,-28},{42,-8}})));
    inner Modelica.Fluid.System system(T_ambient=298.15)
      annotation (Placement(transformation(extent={{60,50},{80,70}})));
  equation
    connect(fixedHeatFlow.port, solution.heatPort) annotation (Line(
        points={{-66,-62},{-58,-62},{-58,-102}},
        color={191,0,0}));
  connect(fixed1.flange, solution.bottom) annotation (Line(
      points={{-18,-84},{2,-84},{2,-102}},
      color={0,127,0}));
    connect(solution.solution, liquidWater.solution) annotation (Line(points={{62,-98},{26,-98},{26,-28}},
                                        color={127,127,0}));
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
    import Chemical;
    extends Modelica.Icons.Example;
    Chemical.Boundaries.Substance Ethanol(
      substanceDefinition=Chemical.Substances.Liquid.Ethanol,
      useFore=false,
      useSolution=true,
      preferMass=true,
      mass_start=(55.508/2)*0.04607) annotation (Placement(transformation(extent={{18,-8},{38,12}})));

    Chemical.Solution solution(useMechanicPorts=true, useThermalPort=true) annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
      annotation (Placement(transformation(extent={{-86,-76},{-66,-56}})));

    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
    Chemical.Boundaries.Substance liquidWater(
      substanceDefinition=Chemical.Substances.Liquid.H2O,
      useFore=false,
      useSolution=true,
      preferMass=true,
      mass_start=1/2) annotation (Placement(transformation(extent={{-50,-8},{-30,12}})));
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

  model ExothermicReaction "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"
    import Chemical;
    extends Modelica.Icons.Example;
    parameter Modelica.Units.SI.MolarEnergy ReactionEnthalpy = -55000;
    Chemical.Solution thermal_isolated_solution(useMechanicPorts = true, ConstantTemperature = false) annotation(
      Placement(transformation(extent = {{-100, -100}, {98, -6}})));
    Chemical.Boundaries.Substance A(useFore = true, useSolution = true, preferMass = false, amountOfSubstance_start = 0.9) annotation(
      Placement(transformation(extent = {{-40, -60}, {-20, -40}})));
    Chemical.Processes.Reaction reaction2_2(process = Chemical.Interfaces.processData(1, ReactionEnthalpy), nS = 1, nP = 1) annotation(
      Placement(transformation(extent = {{-8, -60}, {12, -40}})));
    Chemical.Boundaries.Substance B(useRear = true, useSolution = true, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{28, -60}, {48, -40}})));
    Chemical.Solution solution_at_constant_temperature(useMechanicPorts = true, useThermalPort = true) annotation(
      Placement(transformation(extent = {{-100, 0}, {98, 94}})));
    Chemical.Boundaries.Substance A1(useFore = true, useSolution = true, preferMass = false, amountOfSubstance_start = 0.9) annotation(
      Placement(transformation(extent = {{-40, 40}, {-20, 60}})));
    Chemical.Processes.Reaction reaction2_1(process = Chemical.Interfaces.processData(1, ReactionEnthalpy), nS = 1, nP = 1) annotation(
      Placement(transformation(extent = {{-8, 40}, {12, 60}})));
    Chemical.Boundaries.Substance B1(useRear = true, useSolution = true, preferMass = false, amountOfSubstance_start = 0.1) annotation(
      Placement(transformation(extent = {{20, 40}, {40, 60}})));
    //  Modelica.SIunits.HeatFlowRate q
    //    "Heat flow to environment to reach constant temperature";
    Modelica.Units.SI.Temperature t "Temperature if the solution is ideally thermal isolated from environment";
    Chemical.Boundaries.Substance H2O(substanceDefinition = Chemical.Substances.Liquid.H2O, useSolution = true, mass_start = 1) annotation(
      Placement(transformation(extent = {{20, 4}, {40, 24}})));
    Chemical.Boundaries.Substance H2O1(substanceDefinition = Chemical.Substances.Liquid.H2O, useSolution = true, mass_start = 1) annotation(
      Placement(transformation(extent = {{20, -94}, {40, -74}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1 annotation(
      Placement(transformation(extent = {{-28, 4}, {-8, 24}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed2 annotation(
      Placement(transformation(extent = {{-26, -96}, {-6, -76}})));
    inner Modelica.Fluid.System system(T_ambient = 298.15) annotation(
      Placement(transformation(extent = {{56, 64}, {76, 84}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T = 298.15) annotation(
      Placement(transformation(extent = {{-88, 26}, {-68, 46}})));
  equation
    //  q = fixedTemperature.port.Q_flow;
    t = thermal_isolated_solution.solution.T;
    connect(B.solution, thermal_isolated_solution.solution) annotation(
      Line(points = {{32, -60}, {32, -64}, {58.4, -64}, {58.4, -99.06}}, color = {127, 127, 0}));
    connect(A.solution, thermal_isolated_solution.solution) annotation(
      Line(points = {{-36, -60}, {-36, -64}, {58.4, -64}, {58.4, -99.06}}, color = {127, 127, 0}));
    connect(B1.solution, solution_at_constant_temperature.solution) annotation(
      Line(points = {{24, 40}, {24, 34}, {58.4, 34}, {58.4, 0.94}}, color = {127, 127, 0}));
    connect(A1.solution, solution_at_constant_temperature.solution) annotation(
      Line(points = {{-36, 40}, {-36, 34}, {58.4, 34}, {58.4, 0.94}}, color = {127, 127, 0}));
    connect(solution_at_constant_temperature.solution, H2O.solution) annotation(
      Line(points = {{58.4, 0.94}, {24, 0.94}, {24, 4}}, color = {127, 127, 0}));
    connect(thermal_isolated_solution.solution, H2O1.solution) annotation(
      Line(points = {{58.4, -99.06}, {24, -99.06}, {24, -94}}, color = {127, 127, 0}));
    connect(solution_at_constant_temperature.bottom, fixed1.flange) annotation(
      Line(points = {{-1, -0.94}, {0, -0.94}, {0, 14}, {-18, 14}}, color = {0, 127, 0}));
    connect(thermal_isolated_solution.bottom, fixed2.flange) annotation(
      Line(points = {{-1, -100.94}, {-1, -86}, {-16, -86}}, color = {0, 127, 0}));
    connect(solution_at_constant_temperature.heatPort, fixedTemperature.port) annotation(
      Line(points = {{-60.4, -0.94}, {-60.4, 36}, {-68, 36}}, color = {191, 0, 0}));
    connect(A1.fore, reaction2_1.substrates[1]) annotation(
      Line(points = {{-20, 50}, {-8, 50}}, color = {158, 66, 200}, thickness = 0.5));
    connect(reaction2_1.products[1], B1.rear) annotation(
      Line(points = {{12, 50}, {20, 50}}, color = {158, 66, 200}, thickness = 0.5));
    connect(A.fore, reaction2_2.substrates[1]) annotation(
      Line(points = {{-20, -50}, {-8, -50}}, color = {158, 66, 200}, thickness = 0.5));
    connect(reaction2_2.products[1], B.rear) annotation(
      Line(points = {{12, -50}, {28, -50}}, color = {158, 66, 200}, thickness = 0.5));
    annotation(
      Documentation(revisions = "<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info = "<html>
<p>Demonstration of exotermic reaction with perfect cooling (i.e. connected fixed temperature to the HeatPort) and thermally insulated (HetPort unconnected). See solution_(...).T</p>
</html>"),
      experiment(StopTime = 10, __Dymola_Algorithm = "Dassl"),
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})));
  end ExothermicReaction;

  model HydrogenCombustion "Hydrogen combustion in piston"
    import Chemical;
    extends Modelica.Icons.Example;

    parameter Modelica.Units.SI.Volume V=0.001 "Initial volume";
   // parameter Modelica.SIunits.Pressure p=100000 "Initial pressure";
    parameter Modelica.Units.SI.Temperature T=298.15
      "Initial temperature";

    parameter Modelica.Units.SI.Area A=0.01 "Cross area of cylinder";

    //p*V=n*R*T
    // parameter Modelica.SIunits.AmountOfSubstance n=p*V/(Modelica.Constants.R*T)
    //   "Initial amount of substances in sulution";
    Chemical.Solution idealGas(
      SurfaceArea=A,
      useMechanicPorts=true,
      useThermalPort=true)                                   annotation (Placement(transformation(extent={{-108,-50},{-8,50}})));
                     // AmbientPressure=p)
    //  volume_start=V,
    Chemical.Boundaries.Substance H2_gas(
      substanceDefinition=Chemical.Substances.Gas.H2,
      useFore=true,
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start=26) annotation (Placement(transformation(extent={{-96,-26},{-76,-6}})));
    Chemical.Boundaries.Substance O2_gas(
      substanceDefinition=Chemical.Substances.Gas.O2,
      useFore=true,
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start=13) annotation (Placement(transformation(extent={{-100,10},{-80,30}})));
    Chemical.Boundaries.Substance H2O_gas(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start=1,
      useRear=true) annotation (Placement(transformation(extent={{-36,-8},{-16,12}})));
    Chemical.Processes.Reaction reaction(
      firstProductFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
      firstProduct=Chemical.Substances.Gas.H2O,
      p={2},
      s={1,2},
      redeclare function uLoss = Chemical.Processes.Internal.Kinetics.fastPotentialLoss,
      nS=2,
      nP=1) annotation (Placement(transformation(extent={{-68,-8},{-48,12}})));
    Modelica.Mechanics.Translational.Components.Spring spring(c=1e6) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-58,64})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=200000)
      annotation (Placement(transformation(extent={{-94,-80},{-74,-60}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature coolerTemperature(T=298.15)
      annotation (Placement(transformation(extent={{-18,-80},{-38,-60}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-58,78})));
    Modelica.Mechanics.Translational.Components.Fixed fixed1
      annotation (Placement(transformation(extent={{-68,-66},{-48,-46}})));

    inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=1e-8)
                                               annotation (Placement(transformation(extent={{54,64},{74,84}})));
  equation
  connect(H2_gas.solution, idealGas.solution) annotation (Line(
      points={{-92,-26},{-28,-26},{-28,-49}},
      color={127,127,0}));
  connect(O2_gas.solution, idealGas.solution) annotation (Line(
      points={{-96,10},{-100,10},{-100,-26},{-28,-26},{-28,-49}},
      color={127,127,0}));
  connect(H2O_gas.solution, idealGas.solution) annotation (Line(
      points={{-32,-8},{-32,-28},{-28,-28},{-28,-49}},
      color={127,127,0}));
    connect(idealGas.surfaceFlange, spring.flange_a) annotation (Line(
        points={{-58,50},{-58,54}},
        color={0,127,0}));
    connect(idealGas.heatPort, thermalConductor.port_a) annotation (Line(
        points={{-88,-51},{-88,-56},{-106,-56},{-106,-70},{-94,-70}},
        color={191,0,0}));
    connect(thermalConductor.port_b, coolerTemperature.port) annotation (Line(
        points={{-74,-70},{-38,-70}},
        color={191,0,0}));
    connect(fixed.flange, spring.flange_b) annotation (Line(
        points={{-58,78},{-58,74}},
        color={0,127,0}));
  connect(idealGas.bottom, fixed1.flange) annotation (Line(
      points={{-58,-51},{-58,-56}},
      color={0,127,0}));
    connect(O2_gas.fore, reaction.substrates[1])
      annotation (Line(
        points={{-80,20},{-74,20},{-74,1.75},{-68,1.75}},
        color={158,66,200},
        thickness=0.5));
    connect(H2_gas.fore, reaction.substrates[2])
      annotation (Line(
        points={{-76,-16},{-74,-16},{-74,2.25},{-68,2.25}},
        color={158,66,200},
        thickness=0.5));
    connect(reaction.products[1], H2O_gas.rear) annotation (Line(
        points={{-48,2},{-36,2}},
        color={158,66,200},
        thickness=0.5));
    annotation ( experiment(StopTime=0.00082, __Dymola_Algorithm="Dassl"),
                                         Documentation(info="<html>
<p>The gaseous reaction of hydrogen combustion: </p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p align=\"center\"><b>2 H<sub>2</sub> + O<sub>2</sub> &lt;-&gt; 2 H<sub>2</sub>O</b></p></td>
<td><p>(1)</p></td>
</tr>
</table>
<p><br>This reaction generates a large amount of energy which can be used for mechanical or thermal purposes. </p>
<p>Building this model using the Chemical library components is easy. First, we drag and drop the library class &lsquo;Components.Solution&rsquo; into the diagram of our new model, labeled &lsquo;idealGas&rsquo; in Figure 4. In parameter dialog of this solution we check &ldquo;useThermalPorts&rdquo; and &ldquo;useMechanicsPorts&rdquo; to enable the thermal and mechanical interface. In the same dialog we need to set the area of the piston (e.g., 1 dm<sup>2</sup>), where the pressure provides the force of the green mechanical port of the uppermost side. The next parameter is the ambient external pressure surrounding the system (e.g., 1 bar). All three chemical substances of the reaction (1) can be added by dragging and dropping the library class &lsquo;Components.Substance&rsquo;. Because this model uses gases, the state of matter must be changed to some gas, such as the ideal gas prepared as &lsquo;Interfaces.IdealGas&rsquo;. The substance data must be selected to define the appropriate substances such as &lsquo;Hydrogen_gas&rsquo;, &lsquo;.Oxygen_gas&rsquo; and &lsquo;.Water_gas&rsquo; in package &lsquo;Examples.Substances&rsquo;. In addition, the initial amounts of substances can be prepared for the ideal solution of hydrogen and oxygen gases at a ratio 2:1 to attain the chemical equation above, with the expectation that at the end of the burning process, only water vapor would be presented. Therefore, the initial values of H<sub>2</sub> particles could be set to 26 mmol and of O<sub>2</sub> particles as 13 mmol. All substances must be connected with the &lsquo;idealGas&rsquo; using the blue colored solution port situated on the bottom side of each substance and solution. Then, the chemical reaction is inserted into the diagram of this model as library class &lsquo;Components.Reaction&rsquo;, and it is set to two substrates (nS=2) with stoichiometry s={2,1} and one product with stoichiometry p={2} to represent the reaction (3). The substances are then connected using violet colored substance connectors with appropriate indexes: H<sub>2</sub> to substrates[1], O<sub>2</sub> to substrates[2] and H<sub>2</sub>O to products[1]. At this point, the model is prepared to simulate the conditions of an unconnected heat port and an unconnected mechanical port. This simulation reaches the theoretical ideal of thermally isolated (zero heat flow from/to the solution) and isobaric (zero force generated on piston) conditions. </p>
<p><br><img src=\"modelica://Examples/Resources/Images/Examples/HydrogenBurning.png\"/></p>
<p><font style=\"color: #222222; \">Mueller, M. A., Kim, T. J., Yetter, R. A., &amp; Dryer, F. L. (1999). Flow reactor studies and kinetic modeling of the H2/O2 reaction.&nbsp;<i>International Journal of Chemical Kinetics</i>,&nbsp;<i>31</i>(2), 113-125.</font></p>
<p><br>However, in the real world, there is always some thermal energy flow from the solution, and this cooling process can be connected using the thermal connector of the Modelica Standard Library 3.2.1. For example, the simple thermal conductor of thermal conductance 2W/K at a constant temperature environment of 25&deg;C is represented in the model. The mechanical power of the engine can be connected to the robust mechanical model. However, in our example we selected only a very strong mechanical spring with a spring constant of 10<sup>6</sup> N/m to stop the motion of the piston in order to generate the pressure. This standard spring component is situated above the solution in the model diagram. The results of this experiment are shown in Figure 1. </p>
</html>"),
      Diagram(coordinateSystem(extent={{-120,-100},{120,100}})));
  end HydrogenCombustion;

  model WaterVaporization "Evaporation of water"
    import Chemical;
     extends Modelica.Icons.Example;

    parameter Modelica.Units.SI.Temperature T_start=273.15
      "Initial temperature";

    Chemical.Solution liquid(temperature_start=T_start, useThermalPort=true) annotation (Placement(transformation(extent={{-98,-98},{-6,-8}})));

    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
    Chemical.Solution gas(
      temperature_start=T_start,
      useThermalPort=true)                                   annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                  /*volume_start(
        displayUnit="l") = 0.001, */
    Chemical.Boundaries.Substance H2O_gaseuous(
      useSolution=true,
      preferMass=false,
      useRear=true) annotation (Placement(transformation(extent={{8,50},{28,70}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                         fixedTemperature
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        origin={84,8})));
    Modelica.Blocks.Sources.ContinuousClock clock(offset=1*T_start)
      annotation (Placement(transformation(extent={{56,60},{76,80}})));
    Chemical.Boundaries.Substance otherSubstances(
      substanceDefinition=Chemical.Substances.Gas.O2,
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start=1) annotation (Placement(transformation(extent={{2,28},{22,48}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(
      G=1e6) annotation (Placement(transformation(extent={{44,-12},{64,8}})));
    Chemical.Boundaries.Substance liquidWater(
      substanceDefinition=Chemical.Substances.Liquid.H2O,
      useFore=true,
      useSolution=true,                                                                   mass_start=1)
      annotation (Placement(transformation(extent={{-28,-62},{-48,-42}})));
    inner Modelica.Fluid.System system(p_ambient=100000, T_ambient=298.15)
      annotation (Placement(transformation(extent={{72,-74},{92,-54}})));
    Chemical.Processes.GasSolubility gasVolatility(
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort,
      product=Chemical.Substances.Gas.H2O,
      k_forward=10)
      annotation (Placement(transformation(
          extent={{10,10},{-10,-10}},
          rotation=180,
          origin={-78,26})));
    inner Chemical.DropOfCommons dropOfCommons(L=1) annotation (Placement(transformation(extent={{-92,72},{-72,92}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(G=1e6)
             annotation (Placement(transformation(extent={{44,-40},{64,-20}})));
  equation

    connect(gas.solution, H2O_gaseuous.solution) annotation (Line(
        points={{27.6,6.9},{12,6.9},{12,50}},
        color={127,127,0}));
  connect(fixedTemperature.T, clock.y) annotation (Line(
      points={{96,8},{104,8},{104,70},{77,70}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(gas.solution, otherSubstances.solution) annotation (Line(
      points={{27.6,6.9},{6,6.9},{6,28}},
      color={127,127,0}));
  connect(fixedTemperature.port, thermalConductor.port_b) annotation (Line(
      points={{74,8},{72,8},{72,-2},{64,-2}},
      color={191,0,0}));
  connect(gas.heatPort, thermalConductor.port_a) annotation (Line(
      points={{-27.6,5.1},{-28,5.1},{-28,-2},{44,-2}},
      color={191,0,0}));
    connect(liquid.solution, liquidWater.solution) annotation (Line(points={{-24.4,-97.1},{-24.4,-96.55},{-32,-96.55},{-32,-62}},
                                                           color={127,127,0}));
    connect(liquidWater.fore,gasVolatility. rear) annotation (Line(
        points={{-48,-52},{-88,-52},{-88,26}},
        color={158,66,200},
        thickness=0.5));
    connect(thermalConductor1.port_b, fixedTemperature.port) annotation (Line(points={{64,-30},{68,-30},{68,-2},{72,-2},{72,8},{74,8}}, color={191,0,0}));
    connect(thermalConductor1.port_a, liquid.heatPort) annotation (Line(points={{44,-30},{0,-30},{0,-102},{-79.6,-102},{-79.6,-98.9}}, color={191,0,0}));
    connect(gasVolatility.solution, gas.solution) annotation (Line(points={{-84,16},{-52,16},{-52,6},{27.6,6},{27.6,6.9}},       color={127,127,0}));
    connect(gasVolatility.fore, H2O_gaseuous.rear)
      annotation (Line(
        points={{-68,26},{-52,26},{-52,60},{8,60}},
        color={158,66,200},
        thickness=0.5));
    annotation (
      experiment(
        StopTime=400,
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
    import Chemical;
     extends Modelica.Icons.Example;

    parameter Modelica.Units.SI.Temperature T_start=273.15 - 50
      "Initial temperature";

    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
    Chemical.Solution gas(
      temperature_start=T_start,
      useThermalPort=true,
      BasePressure=600) annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                  /*volume_start(
        displayUnit="l") = 0.001, */
    Chemical.Boundaries.Substance H2O_gaseuous(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start=0.001,
      useRear=true) annotation (Placement(transformation(extent={{-20,56},{0,76}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                         fixedTemperature
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        origin={84,8})));
    Modelica.Blocks.Sources.ContinuousClock clock(offset=1*T_start)
      annotation (Placement(transformation(extent={{62,36},{82,56}})));
    Chemical.Boundaries.Substance otherSubstances(
      substanceDefinition=Chemical.Substances.Gas.O2,
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start=1) annotation (Placement(transformation(extent={{-4,36},{16,56}})));
    Chemical.Solution solid(
      temperature_start=T_start,
      BasePressure=600,
      useThermalPort=true) annotation (Placement(transformation(extent={{-50,-100},{42,-10}})));
    Chemical.Boundaries.Substance H2O_solid(
      substanceDefinition=Chemical.Substances.Solid.H2O_IceIh,
      useFore=true,
      useSolution=true,
      preferMass=true,
      mass_start=1)                   "Solid water" annotation (Placement(transformation(extent={{10,-52},{-10,-32}})));
    Chemical.Processes.GasSolubility gasVolatility(
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort,
      product=Chemical.Substances.Gas.H2O,
      k_forward=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-66,26})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=1e6)
      annotation (Placement(transformation(extent={{48,-8},{68,12}})));
    inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{-88,74},{-68,94}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(G=1e6)
      annotation (Placement(transformation(extent={{50,-32},{70,-12}})));
  equation

    connect(gas.solution, H2O_gaseuous.solution) annotation (Line(
        points={{27.6,6.9},{27.6,0},{-50,0},{-50,56},{-16,56}},
        color={127,127,0}));
  connect(fixedTemperature.T, clock.y) annotation (Line(
      points={{96,8},{98,8},{98,46},{83,46}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(gas.solution, otherSubstances.solution) annotation (Line(
      points={{27.6,6.9},{24,6.9},{24,6},{20,6},{20,36},{0,36}},
      color={127,127,0}));
    connect(solid.solution, H2O_solid.solution) annotation (Line(
        points={{23.6,-99.1},{23.6,-104},{-50,-104},{-50,-52},{6,-52}},
        color={127,127,0}));
    connect(fixedTemperature.port, thermalConductor.port_b) annotation (Line(
        points={{74,8},{72,8},{72,2},{68,2}},
        color={191,0,0}));
    connect(gas.heatPort, thermalConductor.port_a) annotation (Line(
        points={{-27.6,5.1},{-28,5.1},{-28,2},{48,2}},
        color={191,0,0}));
    connect(thermalConductor1.port_a, solid.heatPort) annotation (Line(points={{50,-22},{44,-22},{44,-90},{-31.6,-90},{-31.6,-100.9}}, color={191,0,0}));
    connect(thermalConductor1.port_b, fixedTemperature.port) annotation (Line(points={{70,-22},{72,-22},{72,8},{74,8}}, color={191,0,0}));
    connect(H2O_solid.fore, gasVolatility.rear) annotation (Line(
        points={{-10,-42},{-76,-42},{-76,26}},
        color={158,66,200},
        thickness=0.5));
    connect(gasVolatility.fore, H2O_gaseuous.rear) annotation (Line(
        points={{-56,26},{-56,66},{-20,66}},
        color={200,66,175},
        thickness=0.5));
    connect(gas.solution, gasVolatility.solution) annotation (Line(points={{27.6,6.9},{27.6,12},{-72,12},{-72,16}},       color={127,127,0}));
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
    import Chemical;
     extends Modelica.Icons.Example;

    Chemical.Solution water_solution_25degC(temperature_start=298.15) annotation (Placement(transformation(extent={{-160,-78},{-68,12}})));
                                        //(amountOfSolution_start=52.3)
    Chemical.Solution water_solution_37degC(temperature_start=310.15) annotation (Placement(transformation(extent={{-52,-80},{42,12}})));
                                     //(amountOfSolution_start=39.7)
    Chemical.Processes.GasSolubility CO2_dissolutionP(
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort,
      productFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
      product=Chemical.Substances.Aqueous.CO2,
      k_forward=1) annotation (Placement(transformation(extent={{-138,42},{-118,62}})));
    //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
    Chemical.Boundaries.Substance CO2_25(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start(displayUnit="mmol") = 0.001,
      useRear=true) "Free dissolved CO2 in water at 25 degC" annotation (Placement(transformation(extent={{-130,-28},{-150,-8}})));

    Chemical.Processes.GasSolubility O2_dissolutionP(
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort,
      product=Chemical.Substances.Aqueous.O2,
      k_forward=1) annotation (Placement(transformation(extent={{-98,40},{-78,60}})));

    Chemical.Boundaries.ExternalGas O2_g_25(
      substanceDefinition=Chemical.Substances.Gas.O2,
      useSolution=true,                                                               PartialPressure(displayUnit="mmHg") = 12665.626804425)
      annotation (Placement(transformation(extent={{-114,74},{-94,94}})));
    Chemical.Boundaries.Substance O2_25(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start(displayUnit="mmol") = 0.0001,
      useRear=true) "Free dissolved O2 in water at 25 degC" annotation (Placement(transformation(extent={{-94,-26},{-114,-6}})));

    Chemical.Processes.GasSolubility CO2_dissolutionE(
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort,
      product=Chemical.Substances.Aqueous.CO2,
      k_forward=1) annotation (Placement(transformation(extent={{-24,40},{-4,60}})));

    Chemical.Boundaries.ExternalGas CO2_g_25(
      substanceDefinition=Chemical.Substances.Gas.CO2,
      useSolution=true,                                                                       PartialPressure(displayUnit="mmHg") = 5332.8954966)
      annotation (Placement(transformation(extent={{-154,74},{-134,94}})));

    Chemical.Boundaries.Substance CO2_37(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start(displayUnit="mmol") = 0.001,
      useRear=true) "Free dissolved CO2 in water at 37degC" annotation (Placement(transformation(extent={{-22,-34},{-42,-14}})));

    Chemical.Processes.GasSolubility O2_dissolutionE_NIST(
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort,
      product=Chemical.Substances.Aqueous.O2,
      k_forward=1) annotation (Placement(transformation(extent={{20,42},{40,62}})));
    Chemical.Boundaries.Substance O2_37(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start(displayUnit="mmol") = 0.0001,
      useRear=true) "Free dissolved O2 in water at 37degC" annotation (Placement(transformation(extent={{18,-34},{-2,-14}})));

    Chemical.Boundaries.Substance water_25(
      substanceDefinition=Chemical.Substances.Liquid.H2O,
      useSolution=true,                                                                mass_start=1)
      annotation (Placement(transformation(extent={{-98,-68},{-78,-48}})));
    Chemical.Boundaries.Substance water_37(
      substanceDefinition=Chemical.Substances.Liquid.H2O,
      useSolution=true,                                                                mass_start=1)
      annotation (Placement(transformation(extent={{10,-70},{30,-50}})));
    Chemical.Boundaries.ExternalGas CO2_g_37(
      substanceDefinition=Chemical.Substances.Gas.CO2,
      useSolution=true,                                                                       PartialPressure(displayUnit="mmHg") = 5332.8954966)
      annotation (Placement(transformation(extent={{-44,68},{-24,88}})));
    Chemical.Boundaries.ExternalGas O2_g_37(
      substanceDefinition=Chemical.Substances.Gas.O2,
      useSolution=true,                                                               PartialPressure(displayUnit="mmHg") = 12665.626804425)
      annotation (Placement(transformation(extent={{-6,68},{14,88}})));
    Chemical.Solution water_solution_37degC1(temperature_start=273.15) annotation (Placement(transformation(extent={{66,-80},{160,12}})));
    Chemical.Processes.GasSolubility CO2_dissolutionE1(
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort,
      product=Chemical.Substances.Aqueous.CO2,
      k_forward=1) annotation (Placement(transformation(extent={{92,44},{112,64}})));
    Chemical.Boundaries.Substance CO2_0(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start(displayUnit="mmol") = 0.001,
      useRear=true) "Free dissolved CO2 in water at 0degC" annotation (Placement(transformation(extent={{96,-34},{76,-14}})));

    Chemical.Processes.GasSolubility O2_dissolutionE_NIST1(
      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort,
      product=Chemical.Substances.Aqueous.O2,
      k_forward=1) annotation (Placement(transformation(extent={{134,42},{154,62}})));
    Chemical.Boundaries.Substance O2_0(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start(displayUnit="mmol") = 0.0001,
      useRear=true) "Free dissolved O2 in water at 0degC" annotation (Placement(transformation(extent={{136,-34},{116,-14}})));

    Chemical.Boundaries.Substance water_0(
      substanceDefinition=Chemical.Substances.Liquid.H2O,
      useSolution=true,                                                               mass_start=1)
      annotation (Placement(transformation(extent={{126,-70},{146,-50}})));
    Chemical.Boundaries.ExternalGas CO2_g_0(
      substanceDefinition=Chemical.Substances.Gas.CO2,
      useSolution=true,                                                                      PartialPressure(displayUnit="mmHg") = 5332.8954966)
      annotation (Placement(transformation(extent={{74,68},{94,88}})));
    Chemical.Boundaries.ExternalGas O2_g_0(
      substanceDefinition=Chemical.Substances.Gas.O2,
      useSolution=true,                                                              PartialPressure(displayUnit="mmHg") = 12665.626804425)
      annotation (Placement(transformation(extent={{112,68},{132,88}})));
    inner Modelica.Fluid.System system(p_ambient=100000)
      annotation (Placement(transformation(extent={{-70,-98},{-50,-78}})));
    Real kH_CO2_25, kH_O2_25;
    inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{38,76},{58,96}})));
  equation

    kH_CO2_25 = CO2_25.c / CO2_g_25.substance.x;
    kH_O2_25 = O2_25.c / O2_g_25.substance.x;
  //  kH_CO2_25 = CO2_25.x / CO2_g_25.x;
  //  kH_CO2_25 = CO2_25.x / CO2_g_25.x;

    connect(CO2_25.solution, water_solution_25degC.solution) annotation (Line(
          points={{-134,-28},{-134,-77.1},{-86.4,-77.1}},
                                                        color={127,127,0}));
    connect(O2_25.solution, water_solution_25degC.solution) annotation (Line(
          points={{-98,-26},{-98,-77.1},{-86.4,-77.1}}, color={127,127,0}));
    connect(CO2_37.solution, water_solution_37degC.solution) annotation (Line(
          points={{-26,-34},{-26,-79.08},{23.2,-79.08}},
                                                       color={127,127,0}));
    connect(O2_37.solution, water_solution_37degC.solution) annotation (Line(
          points={{14,-34},{14,-79.08},{23.2,-79.08}}, color={127,127,0}));
    connect(water_25.solution, water_solution_25degC.solution) annotation (Line(
          points={{-94,-68},{-94,-77.1},{-86.4,-77.1}}, color={127,127,0}));
    connect(water_37.solution, water_solution_37degC.solution) annotation (Line(
          points={{14,-70},{14,-79.08},{23.2,-79.08}}, color={127,127,0}));
    connect(CO2_0.solution, water_solution_37degC1.solution) annotation (Line(
          points={{92,-34},{92,-79.08},{141.2,-79.08}},   color={127,127,0}));
    connect(O2_0.solution, water_solution_37degC1.solution) annotation (Line(
          points={{132,-34},{132,-79.08},{141.2,-79.08}}, color={127,127,0}));
    connect(water_0.solution, water_solution_37degC1.solution) annotation (Line(
          points={{130,-70},{130,-79.08},{141.2,-79.08}}, color={127,127,0}));
    connect(CO2_g_25.fore, CO2_dissolutionP.rear) annotation (Line(
        points={{-134,84},{-138,84},{-138,52}},
        color={158,66,200},
        thickness=0.5));
    connect(O2_g_25.fore, O2_dissolutionP.rear) annotation (Line(
        points={{-94,84},{-98,84},{-98,50}},
        color={158,66,200},
        thickness=0.5));
    connect(CO2_g_37.fore, CO2_dissolutionE.rear) annotation (Line(
        points={{-24,78},{-24,50}},
        color={158,66,200},
        thickness=0.5));
    connect(O2_g_37.fore, O2_dissolutionE_NIST.rear) annotation (Line(
        points={{14,78},{20,78},{20,52}},
        color={158,66,200},
        thickness=0.5));
    connect(CO2_g_0.fore, CO2_dissolutionE1.rear) annotation (Line(
        points={{94,78},{92,78},{92,54}},
        color={158,66,200},
        thickness=0.5));
    connect(O2_g_0.fore, O2_dissolutionE_NIST1.rear) annotation (Line(
        points={{132,78},{132,52},{134,52}},
        color={158,66,200},
        thickness=0.5));
    connect(CO2_g_25.solution, water_solution_25degC.solution)
      annotation (Line(points={{-150,74},{-150,70},{-156,70},{-156,-78},{-122,-78},{-122,-77.1},{-86.4,-77.1}}, color={127,127,0}));
    connect(O2_g_25.solution, water_solution_25degC.solution)
      annotation (Line(points={{-110,74},{-110,68},{-156,68},{-156,-78},{-122,-78},{-122,-77.1},{-86.4,-77.1}}, color={127,127,0}));
    connect(CO2_g_37.solution, water_solution_37degC.solution) annotation (Line(points={{-40,68},{-46,68},{-46,-79.08},{23.2,-79.08}}, color={127,127,0}));
    connect(O2_g_37.solution, water_solution_37degC.solution) annotation (Line(points={{-2,68},{-46,68},{-46,-79.08},{23.2,-79.08}}, color={127,127,0}));
    connect(CO2_g_0.solution, water_solution_37degC1.solution) annotation (Line(points={{78,68},{70,68},{70,-79.08},{141.2,-79.08}}, color={127,127,0}));
    connect(O2_g_0.solution, water_solution_37degC1.solution) annotation (Line(points={{116,68},{70,68},{70,-79.08},{141.2,-79.08}}, color={127,127,0}));
    connect(CO2_dissolutionP.fore, CO2_25.rear) annotation (Line(
        points={{-118,52},{-118,-18},{-130,-18}},
        color={200,66,175},
        thickness=0.5));
    connect(O2_dissolutionP.fore, O2_25.rear) annotation (Line(
        points={{-78,50},{-78,-16},{-94,-16}},
        color={200,66,175},
        thickness=0.5));
    connect(CO2_dissolutionE.fore, CO2_37.rear) annotation (Line(
        points={{-4,50},{-4,-24},{-22,-24}},
        color={200,66,175},
        thickness=0.5));
    connect(O2_dissolutionE_NIST.fore, O2_37.rear) annotation (Line(
        points={{40,52},{40,-24},{18,-24}},
        color={200,66,175},
        thickness=0.5));
    connect(CO2_dissolutionE1.fore, CO2_0.rear)
      annotation (Line(
        points={{112,54},{112,-26},{96,-26},{96,-24}},
        color={200,66,175},
        thickness=0.5));
    connect(O2_dissolutionE_NIST1.fore, O2_0.rear)
      annotation (Line(
        points={{154,52},{154,-26},{136,-26},{136,-24}},
        color={200,66,175},
        thickness=0.5));
    connect(CO2_dissolutionP.solution, water_solution_25degC.solution)
      annotation (Line(points={{-134,42},{-134,34},{-76,34},{-76,-78},{-86.4,-78},{-86.4,-77.1}},       color={127,127,0}));
    connect(water_solution_25degC.solution, O2_dissolutionP.solution)
      annotation (Line(points={{-86.4,-77.1},{-86.4,-78},{-76,-78},{-76,40},{-94,40}},       color={127,127,0}));
    connect(water_solution_37degC.solution, CO2_dissolutionE.solution)
      annotation (Line(points={{23.2,-79.08},{23.2,-80},{40,-80},{40,36},{-2,36},{-2,40},{-20,40}},      color={127,127,0}));
    connect(water_solution_37degC.solution, O2_dissolutionE_NIST.solution)
      annotation (Line(points={{23.2,-79.08},{23.2,-80},{40,-80},{40,42},{24,42}},       color={127,127,0}));
    connect(CO2_dissolutionE1.solution, water_solution_37degC1.solution)
      annotation (Line(points={{96,44},{96,34},{168,34},{168,-79.08},{141.2,-79.08}},         color={127,127,0}));
    connect(O2_dissolutionE_NIST1.solution, water_solution_37degC1.solution)
      annotation (Line(points={{138,42},{168,42},{168,-79.08},{141.2,-79.08}},       color={127,127,0}));
    annotation (
      experiment(StopTime=1, __Dymola_Algorithm="Dassl"),
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

  model EnzymeKinetics "Basic enzyme kinetics"
    extends Modelica.Icons.Example;
    Chemical.Solution solution annotation(
      Placement(transformation(extent = {{-100, -100}, {100, 100}})));
    //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
    Chemical.Boundaries.Substance P(useRear = true, useSolution = true, preferMass = false, amountOfSubstance_start = 1e-8) annotation(
      Placement(transformation(extent = {{72, -12}, {92, 8}})));
    Chemical.Boundaries.Substance S(useFore = true, useSolution = true, preferMass = false, amountOfSubstance_start = 100) annotation(
      Placement(transformation(extent = {{-92, -14}, {-72, 6}})));
    parameter Modelica.Units.SI.AmountOfSubstance tE = 1 "Total amount of enzyme";
    parameter Real k_cat(unit = "mol/s", displayUnit = "mol/min") = 1 "Forward rate of second reaction";
    constant Modelica.Units.SI.Concentration Km = 0.1 "Michaelis constant = substrate concentration at rate of half Vmax";
    parameter Modelica.Units.SI.MolarFlowRate Vmax = 1e-5*k_cat "Maximal molar flow";
    Chemical.Boundaries.Substance ES(useRear = true, useFore = true, useSolution = true, preferMass = false,
      amountOfSubstance_start=tE/1e2)                                                                                                         annotation(
      Placement(transformation(extent = {{-8, -10}, {12, 10}})));
    Chemical.Boundaries.Substance E(useRear = true, useFore = true, useSolution = true, preferMass = false,
      amountOfSubstance_start=tE*(1 - 1e-2))                                                                                                     annotation(
      Placement(transformation(extent = {{12, 36}, {-8, 56}})));
    Processes.Reaction chemicalReaction(process = Chemical.Interfaces.processData(2/Km), k_forward = 1, nS = 2, nP = 1) annotation(
      Placement(transformation(extent = {{-42, -8}, {-22, 12}})));
    Processes.ForwardReaction chemicalReaction1(k_forward = k_cat, nS = 1, nP = 2) annotation(
      Placement(transformation(extent = {{26, -8}, {46, 12}})));
    Chemical.Boundaries.Substance liquidWater(useSolution = true, substanceDefinition = Chemical.Substances.Liquid.H2O, preferMass = true, mass_start = 1) annotation(
      Placement(transformation(extent = {{42, -80}, {62, -60}})));
    inner DropOfCommons dropOfCommons(n_flow_reg=1e-8)
                                      annotation(
      Placement(transformation(extent = {{68, 70}, {88, 90}})));
  equation
    //Michaelis-Menton: v=((E.q_out.conc + ES.q_out.conc)*k_cat)*S.concentration/(Km+S.concentration);
    connect(E.solution, solution.solution) annotation(
      Line(points = {{8, 36}, {-8, 36}, {-8, -98}, {60, -98}}, color = {127, 127, 0}));
    connect(ES.solution, solution.solution) annotation(
      Line(points = {{-4, -10}, {-4, -98}, {60, -98}}, color = {127, 127, 0}));
    connect(S.solution, solution.solution) annotation(
      Line(points = {{-88, -14}, {-88, -56}, {-8, -56}, {-8, -98}, {60, -98}}, color = {127, 127, 0}));
    connect(P.solution, solution.solution) annotation(
      Line(points = {{76, -12}, {76, -98}, {60, -98}}, color = {127, 127, 0}));
    connect(liquidWater.solution, solution.solution) annotation(
      Line(points = {{46, -80}, {46, -98}, {60, -98}}, color = {127, 127, 0}));
    connect(chemicalReaction.products[1], ES.rear) annotation(
      Line(points = {{-22, 2}, {-16, 2}, {-16, 0}, {-8, 0}}, color = {158, 66, 200}, thickness = 0.5));
    connect(ES.fore, chemicalReaction1.substrates[1]) annotation(
      Line(points = {{12, 0}, {18, 0}, {18, 2}, {26, 2}}, color = {158, 66, 200}, thickness = 0.5));
    connect(S.fore, chemicalReaction.substrates[1]) annotation(
      Line(points = {{-72, -4}, {-52, -4}, {-52, -2}, {-42, -2}, {-42, 1.75}}, color = {158, 66, 200}, thickness = 0.5));
    connect(E.fore, chemicalReaction.substrates[2]) annotation(
      Line(points = {{-8, 46}, {-52, 46}, {-52, 2.25}, {-42, 2.25}}, color = {158, 66, 200}, thickness = 0.5));
    connect(chemicalReaction1.products[1], P.rear) annotation(
      Line(points = {{46, 1.75}, {58, 1.75}, {58, -2}, {72, -2}}, color = {158, 66, 200}, thickness = 0.5));
    connect(E.rear, chemicalReaction1.products[2]) annotation(
      Line(points = {{12, 46}, {48, 46}, {48, 48}, {56, 48}, {56, 2.25}, {46, 2.25}}, color = {158, 66, 200}, thickness = 0.5));
    annotation(
      Documentation(revisions = "<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info = "<html>
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
      experiment(StopTime=70000, __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput);
  end EnzymeKinetics;

  model WaterElectrolysis "Water electrolysis"
    import Chemical;

    extends Modelica.Icons.Example;
    Chemical.Boundaries.Substance O2_gas(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start=0.001,
      useRear=true) annotation (Placement(transformation(extent={{-12,-6},{8,14}})));

    Chemical.Boundaries.Substance H2_gas(
      useSolution=true,
      preferMass=false,
      amountOfSubstance_start=0.001,
      useRear=true) annotation (Placement(transformation(extent={{38,-6},{18,14}})));
    Chemical.Processes.Reaction reaction(
      s={2,4},
      p={4,2,1},
      firstProductFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
      firstProduct=Chemical.Substances.Solid.e,
      nextProducts={Chemical.Substances.Aqueous.H2,Chemical.Substances.Aqueous.O2},
      nP=3,
      nS=2) annotation (Placement(transformation(
          extent={{11,11},{-11,-11}},
          rotation=180,
          origin={-31,-29})));
    Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{58,-78},{92,30}})));
    Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-90,-80},{-56,28}})));
    Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
      annotation (Placement(transformation(extent={{-42,70},{-22,90}})));
    Chemical.Boundaries.ElectronTransfer electrone(useRear=true, useFore=false) annotation (Placement(transformation(extent={{66,-38},{86,-18}})));
    Chemical.Boundaries.ElectronTransfer electrone1(useFore=true) annotation (Placement(transformation(extent={{-84,-34},{-64,-14}})));
  Modelica.Electrical.Analog.Basic.Resistor resistor(R=1)
    annotation (Placement(transformation(extent={{-36,38},{-16,58}})));
  Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
    annotation (Placement(transformation(extent={{-66,38},{-46,58}})));
    Chemical.Solution air                                                                 annotation (Placement(transformation(extent={{-40,-16},{50,26}})));
    Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
      annotation (Placement(transformation(extent={{18,38},{-2,58}})));
    Chemical.Boundaries.Substance liquidWater(
      substanceDefinition=Chemical.Substances.Liquid.H2O,
      useFore=true,
      useSolution=false,                                                                  mass_start=1)
      annotation (Placement(transformation(extent={{0,-72},{-20,-52}})));
    Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{36,26},{56,46}})));
  equation
  connect(electrone1.pin,voltageSensor. p) annotation (Line(
      points={{-74,-14.2},{-92,-14.2},{-92,48},{-74,48},{-74,80},{-42,80}},
      color={0,0,255}));
  connect(electrone.pin,voltageSensor. n) annotation (Line(
      points={{76,-18.2},{76,48},{60,48},{60,80},{-22,80}},
      color={0,0,255}));
  connect(electrone.solution,anode. solution) annotation (Line(
      points={{70,-38},{70,-68},{84,-68},{84,-76.92},{85.2,-76.92}},
      color={127,127,0}));
  connect(electrone1.pin,currentSensor. p) annotation (Line(
      points={{-74,-14.2},{-92,-14.2},{-92,48},{-66,48}},
      color={0,0,255}));
  connect(currentSensor.n,resistor. p) annotation (Line(
      points={{-46,48},{-36,48}},
      color={0,0,255}));
  connect(electrone1.solution,cathode. solution) annotation (Line(
      points={{-80,-34},{-80,-70},{-64,-70},{-64,-78.92},{-62.8,-78.92}},
      color={127,127,0}));
    connect(O2_gas.solution, air.solution) annotation (Line(points={{-8,-6},{-8,-20},{32,-20},{32,-15.58}},
                                           color={127,127,0}));
    connect(H2_gas.solution, air.solution) annotation (Line(points={{34,-6},{34,-15.58},{32,-15.58}},
                          color={127,127,0}));
    connect(constantVoltage.p, voltageSensor.n) annotation (Line(points={{18,48},{
            60,48},{60,80},{-22,80}}, color={0,0,255}));
    connect(resistor.n, constantVoltage.n)
      annotation (Line(points={{-16,48},{-2,48}}, color={0,0,255}));
    connect(ground.p, constantVoltage.p) annotation (Line(points={{46,46},{46,48},{18,48}}, color={0,0,255}));
    connect(electrone.rear, reaction.products[1])
      annotation (Line(
        points={{66,-28},{62,-28},{62,-30},{-18,-30},{-18,-29.3667},{-20,-29.3667}},
        color={158,66,200},
        thickness=0.5));
    connect(H2_gas.rear, reaction.products[2]) annotation (Line(
        points={{38,4},{54,4},{54,-29},{-20,-29}},
        color={158,66,200},
        thickness=0.5));
    connect(O2_gas.rear, reaction.products[3])
      annotation (Line(
        points={{-12,4},{-12,-28.6333},{-20,-28.6333}},
        color={158,66,200},
        thickness=0.5));
    connect(liquidWater.fore, reaction.substrates[1])
      annotation (Line(
        points={{-20,-62},{-52,-62},{-52,-29.275},{-42,-29.275}},
        color={158,66,200},
        thickness=0.5));
    connect(electrone1.fore, reaction.substrates[2])
      annotation (Line(
        points={{-64,-24},{-48,-24},{-48,-28.725},{-42,-28.725}},
        color={158,66,200},
        thickness=0.5));
    annotation ( experiment(StopTime=1), Documentation(info="<html>
<p>The water ecectrolysis: </p>
<p><b>2 H<sub>2</sub>O +&nbsp;&nbsp;4 e<sup>-</sup><sub>(catode)</sub>&nbsp;&lt;-&gt;  2 H<sub>2</sub> + O<sub>2</sub>&nbsp;+&nbsp;&nbsp;4 e<sup>-</sup><sub>(anode)</sub>&nbsp;</b></p>
</html>"));
  end WaterElectrolysis;

  model SimpleReactionsWithJunction
    extends Modelica.Icons.Example;
    Chemical.Processes.Reaction r( nP = 1,         process = Chemical.Interfaces.processData(2),
      nS=1)                                                                                      annotation(
      Placement(transformation(extent={{-6,12},{14,32}})));
    Chemical.Boundaries.Substance A(useFore = true) annotation(
      Placement(transformation(extent={{-68,-10},{-48,10}})));
    Chemical.Boundaries.Substance B(useRear = true) annotation(
      Placement(transformation(extent = {{30, 14}, {50, 34}})));
    Boundaries.Substance          B1(useRear=true)  annotation(
      Placement(transformation(extent={{34,-32},{54,-12}})));
    Processes.Reaction          r1(
      nP=1,
      nS=1,
      process=Chemical.Interfaces.processData(2))                                                annotation(
      Placement(transformation(extent={{-2,-32},{18,-12}})));
    Processes.Reaction r_ref(
      nP=1,
      nS=1,
      process=Chemical.Interfaces.processData(2)) annotation (Placement(transformation(extent={{-2,64},{18,84}})));
    Boundaries.Substance A_ref(useFore=true) annotation (Placement(transformation(extent={{-54,64},{-34,84}})));
    Boundaries.Substance B_ref(useRear=true) annotation (Placement(transformation(extent={{34,66},{54,86}})));
    Topology.JunctionRFF junctionRFF annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
  equation
    connect(r.products[1], B.rear) annotation(
      Line(points={{14,22},{22,22},{22,24},{30,24}},
                                          color = {158, 66, 200}, thickness = 0.5));
    connect(r1.products[1], B1.rear) annotation (Line(
        points={{18,-22},{34,-22}},
        color={158,66,200},
        thickness=0.5));
    connect(A_ref.fore, r_ref.substrates[1]) annotation (Line(
        points={{-34,74},{-2,74}},
        color={158,66,200},
        thickness=0.5));
    connect(r_ref.products[1], B_ref.rear) annotation (Line(
        points={{18,74},{26,74},{26,76},{34,76}},
        color={158,66,200},
        thickness=0.5));
    connect(A.fore, junctionRFF.rear) annotation (Line(
        points={{-48,0},{-30,0}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionRFF.foreA, r.substrates[1]) annotation (Line(
        points={{-20,10},{-20,22},{-6,22}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionRFF.foreB, r1.substrates[1]) annotation (Line(
        points={{-20,-10},{-20,-22},{-2,-22}},
        color={158,66,200},
        thickness=0.5));
    annotation(
      Documentation(revisions = "<html>
<p><i>2025</i></p>
<p>Marek Matejak</p>
</html>", info = "<html>
        <p>Simple reaction demonstrating equilibration between substance A and substance B, in constant solution. Observe the molar concentration (A.c) and molar fraction.</p>
</html>"),
      experiment(StopTime = 10));
  end SimpleReactionsWithJunction;

  annotation(
    Documentation(info = "<html>
<u>Tests for top level components of the undirected chemical simulation package.</u>
</html>"));
end Examples;
