within Chemical;
package Undirected
  package Interfaces
    connector Fore "Undirected connector outputting the forward state"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
      flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

      output Chemical.Interfaces.SubstanceState state_forwards "State of substance in forwards direction";
      input Chemical.Interfaces.SubstanceState state_rearwards "State of substance in rearwards direction";

      output Chemical.Interfaces.SolutionState solution "State of solution";
      output stateOfMatter.SubstanceData definition "Definition of substance";
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Ellipse(
              extent={{-80,80},{80,-80}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Ellipse(
              extent={{-40,40},{40,-40}},
              lineColor={194,138,221},
              lineThickness=0.5,
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid)}),
                                Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>One of the two undirected connectors. The state information flows in both directions, forward and backward. </u>
<u>At positive molarflow, Fore is a output. </u>
</html>"));

    end Fore;

    connector Rear "Undirected connector outputting the rearward state"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
      flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

      input Chemical.Interfaces.SubstanceState state_forwards "State of substance in forwards direction";
      output Chemical.Interfaces.SubstanceState state_rearwards "State of substance in rearwards direction";

      input Chemical.Interfaces.SolutionState solution "State of solution";
      input stateOfMatter.SubstanceData definition "Definition of substance";

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Ellipse(
              extent={{-80,80},{80,-80}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid), Ellipse(
              extent={{-40,40},{40,-40}},
              lineColor={194,138,221},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
                                Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>One of the two undirected connectors. The state information flows in both directions, forward and backward. </u>
<u>At positive molarflow, Rear is a input.</u>
</html>"));
    end Rear;

    partial model SISOBiFlow "Base Model with basic flow eqautions for SISO"

      import Chemical.Utilities.Types.InitializationMethods;

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the flow" annotation (Dialog(tab="Advanced"));
      parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
        annotation(Dialog(tab="Advanced"));
      parameter InitializationMethods initN_flow = Chemical.Utilities.Types.InitializationMethods.none "Initialization method for n_flow"
        annotation(Dialog(tab= "Initialization"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_0=0 "Initial value for n_flow"
        annotation (Dialog(tab="Initialization", enable=(initM_flow == InitializationMethods.state)));
      parameter Chemical.Utilities.Units.MolarFlowAcceleration n_acceleration_0=0 "Initial value for der(n_flow)"
        annotation (Dialog(tab="Initialization", enable=(initM_flow == InitializationMethods.derivative)));

      Chemical.Undirected.Interfaces.Fore fore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));
      Chemical.Undirected.Interfaces.Rear rear(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

      Modelica.Units.SI.MolarFlowRate n_flow(stateSelect=n_flowStateSelect)=rear.n_flow;

      // delta potential computation
      Modelica.Units.SI.ChemicalPotential du_fore;
                           // delta = fore - rear
      Modelica.Units.SI.ChemicalPotential du_rear;
                           // delta = rear - fore

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      // input state quantities
      Modelica.Units.SI.ChemicalPotential u_rear_in=rear.state_forwards.u "Chemical potential of substance entering";
      Modelica.Units.SI.ChemicalPotential u_fore_in=fore.state_rearwards.u "Chemical potential of substance entering";
      Modelica.Units.SI.MolarEnthalpy h_rear_in=rear.state_forwards.h "Enthalpy of substance enetering";
      Modelica.Units.SI.MolarEnthalpy h_fore_in=fore.state_rearwards.h "Enthalpy of substance enetering";

      //outlet state quantities
      Modelica.Units.SI.ChemicalPotential u_rear_out "Chemical potential of substance exiting";
      Modelica.Units.SI.ChemicalPotential u_fore_out "Chemical potential of substance exiting";
      Modelica.Units.SI.MolarEnthalpy h_rear_out "Enthalpy of substance exiting";
      Modelica.Units.SI.MolarEnthalpy h_fore_out "Enthalpy of substance exiting";

    initial equation
      if initN_flow == InitializationMethods.state then
        n_flow = n_flow_0;
      elseif initN_flow == InitializationMethods.derivative then
        der(n_flow) = n_acceleration_0;
      elseif initN_flow == InitializationMethods.steadyState then
        der(n_flow) = 0;
      end if;

    equation

      fore.n_flow + rear.n_flow = 0;
      fore.r = rear.r - der(rear.n_flow) * L;

      u_fore_out = u_rear_in + du_fore;
      u_rear_out = u_fore_in + du_rear;


      rear.state_rearwards.u = u_rear_out;
      rear.state_rearwards.h = h_rear_out;
      fore.state_forwards.u = u_fore_out;
      fore.state_forwards.u = h_fore_out;

      annotation (Documentation(info="<html>
<u>Interface class for all components with one fore and one rear port and a massflow without a mass storage between.</u>
<u>This class already implements the equations that are common for such components, namly the conservation of mass, the intertance equation, as well as the clipping of u_out to u_min. </u>
<u>If u_out should be lower the u_min, the remaining potential drop is added on the difference in inertial potential r, basically accelerating or decelerating the massflow. </u>
<u>The component offers different initialization methods for the massflow, as well as several parameters used in the equations above. </u>
<u>The clipping of the massflow can be turned off (this should be done by the modeler as a final modificator while extending to hide this option from the enduser).</u>
</html>"));
    end SISOBiFlow;
  end Interfaces;

  package Boundaries "Boundary models for undirected chemical simulation"
    extends Modelica.Icons.SourcesPackage;

    model Substance "Substance in solution"
      extends Icons.Substance;
      extends Internal.PartialSubstanceInSolution(
        useSolution=false,
        useFore=false,
        useRear=false);

      import Chemical.Utilities.Types.InitializationMethods;

      parameter InitializationMethods initAmount = Chemical.Utilities.Types.InitializationMethods.state "Initialization method for amount of substance"
       annotation (HideResult=not useRear, Dialog(enable=useRear));
       // annotation(Dialog(tab= "Initialization", group="Molar flow"));


      Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
        "Molar concentration of particles";

       parameter stateOfMatter.SubstanceDataParameters substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true, Dialog(enable=not useRear));

      parameter Boolean use_mass_start = true "use mass_start, otherwise amountOfSubstance_start"
        annotation (Evaluate=true, choices(checkBox=true), Dialog(group="Initialization"));

    parameter Modelica.Units.SI.Mass mass_start=1
      "Initial mass of the substance"
      annotation (HideResult=not use_mass_start, Dialog(group="Initialization", enable=use_mass_start));

    parameter Modelica.Units.SI.AmountOfSubstance amountOfSubstance_start=1
      "Initial amount of substance base molecules"
        annotation (HideResult=use_mass_start, Dialog(group="Initialization", enable=not use_mass_start));

      Modelica.Units.SI.Mass mass=amountOfBaseMolecules*
          molarMassOfBaseMolecule "Mass";

      parameter Boolean calculateClusteringHeat = true "Only for self clustering substances"
          annotation(Evaluate=true, choices(checkBox=true), Dialog(tab = "Clustering", enable = stateOfMatter.selfClustering(substanceData)));

    protected
      //parameter
      Modelica.Units.SI.Mass m_start=if use_mass_start then mass_start else
         amountOfSubstance_start*molarMassOfBaseMolecule;

      //parameter
      Modelica.Units.SI.MolarMass molarMassOfBaseMolecule = stateOfMatter.molarMassOfBaseMolecule(substanceDataVar);

      Modelica.Units.SI.AmountOfSubstance amountOfBaseMolecules(start=amountOfSubstance_start)
        "Amount of base molecules inside all clusters in compartment";

      Modelica.Units.SI.AmountOfSubstance amountOfFreeMolecule(start=
           m_start*stateOfMatter.specificAmountOfFreeBaseMolecule(
                                       substanceDataVar,
                                       T=system.T_ambient,
                                       p=system.p_ambient))
        "Amount of free molecules not included inside any clusters in compartment";

      Modelica.Units.SI.AmountOfSubstance amountOfParticles(start=
           m_start*stateOfMatter.specificAmountOfParticles(
                                       substanceDataVar,
                                       T=system.T_ambient,
                                       p=system.p_ambient))
        "Amount of particles/clusters in compartment";

      Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
          Modelica.Constants.R*solutionState.T))
        "Dissociation constant of hydrogen bond between base molecules";

      Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
          stateOfMatter.selfClusteringBondEnthalpy(substanceDataVar)
        - solutionState.T * stateOfMatter.selfClusteringBondEntropy(substanceDataVar)
        "Gibbs energy of hydrogen bond between H2O molecules";

      Modelica.Units.SI.AmountOfSubstance amountOfBonds
        "Amount of hydrogen bonds between molecules in compartment";

      Real logn(stateSelect=StateSelect.prefer, start=log(m_start/molarMassOfBaseMolecule), min=0)
      "Natural logarithm of the amount of base molecules in solution";

      parameter Modelica.Units.SI.MolarFlowRate change_start=1
      "Initial change of substance base molecules"
        annotation ( Dialog(group="Initialization", enable=(initAmount == InitializationMethods.derivative)));


      parameter Boolean EnthalpyNotUsed=false annotation (
        Evaluate=true,
        HideResult=true,
        choices(checkBox=true),
        Dialog(tab="Advanced", group="Performance"));

    initial equation

      if not useRear then
        amountOfBaseMolecules = m_start/molarMassOfBaseMolecule;
      elseif initAmount == InitializationMethods.steadyState then
        r_fore_intern=0;
      elseif initAmount == InitializationMethods.state then
        amountOfBaseMolecules = amountOfSubstance_start;
      elseif initAmount == InitializationMethods.derivative then
        n_flow = change_start;
      end if;


    equation
      if not useRear then
       substanceDataVar = substanceData;
      end if;
     //n_flow = n_flow_out;

      if stateOfMatter.selfClustering(substanceDataVar) then

        //Liquid cluster theory - equilibrium:
        //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules
        //amountOfParticles/solutionState.n = x/(1-K*x);                //sum(x[i])
        //amountOfBaseMolecules/solutionState.n = x/((1-K*x)^2);            //sum(i*x[i])
        //amountOfHydrogenBonds/solutionState.n = x*x*K/((1-K*x)^2);   //sum((i-1)*x[i])

        amountOfParticles*(1 - SelfClustering_K*substance.x) = amountOfFreeMolecule;

        //Calculation of "abs(amountOfBaseMolecules*(1 - SelfClustering_K*x)) = amountOfParticles":
        substance.x = ((2*SelfClustering_K+solutionState.n/amountOfBaseMolecules) - sqrt((4*SelfClustering_K*solutionState.n/amountOfBaseMolecules)+(solutionState.n/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        amountOfBonds = amountOfBaseMolecules*substance.x*SelfClustering_K;

        //TODO: may be the volume of the same number of free water molecules is different as volume of the same number of water molecules in cluster ..
        //TODO: more precise calculation of other properties

       //der(enthalpy) = solutionState.dH + n_flow*actualStream(port_a.h_outflow);
       //enthalpy = molarEnthalpy*amountOfBaseMolecules + amountOfAdditionalBonds*bondEnthalpy;
        dH =if (EnthalpyNotUsed) then 0 else der(substance.h)*
          amountOfBaseMolecules +
          n_flow*substance.h - h_flow +
          (if (calculateClusteringHeat) then stateOfMatter.selfClusteringBondEnthalpy(
          substanceDataVar)*der(amountOfBonds) else 0)
                        "heat transfer from other substances in solution [J/s]";

        Gj =amountOfBaseMolecules*substance.u + amountOfBonds*SelfClustering_dG
                        "Gibbs energy of the substance";

      else

        amountOfParticles = amountOfFreeMolecule;
        amountOfBaseMolecules = amountOfFreeMolecule;
        amountOfBonds = 0;

        //der(enthalpy) = solutionState.dH + n_flow*actualStream(port_a.h_outflow);
        //enthalpy = molarEnthalpy*amountOfBaseMolecules;
        dH =
          if (EnthalpyNotUsed) then  0
          else    der(substance.h)*amountOfBaseMolecules +
                  n_flow*substance.h - h_flow
                  "heat transfer from other substances in solution [J/s]";

        Gj = amountOfBaseMolecules*substance.u "Gibbs energy of the substance [J]";

      end if;

      //The main accumulation equation is "der(amountOfBaseMolecules)=n_flow"
      // However, the numerical solvers can handle it in form of log(n) much better. :-)
      der(logn) = (n_flow/amountOfBaseMolecules) "accumulation of amountOfBaseMolecules=exp(logn) [mol]";
      //der(amountOfBaseMolecules) = n_flow;
      amountOfBaseMolecules = exp(logn);

      substance.x = amountOfFreeMolecule/solutionState.n "mole fraction [mol/mol]";

      c = amountOfParticles/solutionState.V "concentration [mol/m3]";



      //solution flows
      i = Modelica.Constants.F*substance.z*n_flow +
          Modelica.Constants.F*der(substance.z)*amountOfBaseMolecules "change of sunstance charge [A]";
      dV = substance.Vm*n_flow + der(substance.Vm)*amountOfBaseMolecules "change of substance volume [m3/s]";

      //extensive properties
      nj = amountOfParticles;
      mj = amountOfBaseMolecules*molarMassOfBaseMolecule;
      Vj = amountOfBaseMolecules*substance.Vm;
      Qj = Modelica.Constants.F*amountOfBaseMolecules*substance.z;
      Ij = (1/2)*(amountOfBaseMolecules*substance.z^2);

         annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={Text(
              extent={{-84,22},{92,64}},
              lineColor={128,0,255},
              textString="%name")}), Documentation(revisions="<html>
<p>2009-2025 by Marek Matejak, Ph.D. </p>
</html>",     info="<html>
<h4>n = x &middot; n(solution) = &int; MolarFlow</h4>
<p>where n is amount of the substance and x is mole fraction.</p>
<p>The main class from &ldquo;Chemical&rdquo; package is called &quot;Substance&quot;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &quot;n&quot; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &quot;n(solution)&quot; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solutionState. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
<p><br>The recalculation between mole fraction, molarity and molality can be written as follows:</p>
<p>x = n/n(solution) = b * m(solvent)/n(solution) = c * V(solution)/n(solution)</p>
<p>where m(solvent) is mass of solvent, V(solution) is volume of solution, b=n/m(solvent) is molality of the substance, c=n/V(solution) is molarity of the substance.</p>
<p>If the amount of solution is selected to the number of total solution moles per one kilogram of solvent then the values of x will be the same as molality.</p>
<p>If the amount of solution is selected to the number of total solution moles in one liter of solution then the values of x will be the same as molarity.</p>
<p><br><br>Definition of electro-chemical potential:</p>
<h4>u = u&deg; + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg; = DfG = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. relative eletric potential of the solution</p>
<p>z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
<p>R .. gas constant</p>
<p>F .. Faraday constant</p>
<p>gamma .. activity coefficient</p>
<p>u&deg; .. chemical potential of pure substance</p>
<p>DfG .. free Gibbs energy of formation of the substance</p>
<p>DfH .. free enthalpy of formation of the substance</p>
<p>DfS .. free entropy of formation of the substance </p>
<p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
</html>"));
    end Substance;

    model ElectronTransfer "Electron transfer from the solution to electric circuit"
      extends Icons.ElectronTransfer;
      extends Internal.PartialSubstanceInSolution(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible,useSolution=false,
        useFore=false,
        useRear=false);

      Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
          Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
              extent={{-10,88},{10,108}})));

      parameter Modelica.Units.SI.ChemicalPotential u_0=0 "Initial electro-chemical potential"
         annotation (HideResult=useRear, Dialog(group="Initialization", enable=not useRear));


    initial equation
      if not useRear then
        substance.u = u_0;
      else
        substance.u = state_in_rear.u;
      end if;

    equation
      if not useRear then
       substanceDataVar = Chemical.Substances.Electrone_solid();
      end if;


      substance.x = 1;

      //electric
      pin.v = substance.v;
      pin.i + substance.z*Modelica.Constants.F*n_flow + i = 0;

      //none solution changes
      dH = 0;
      dV = 0;
      nj=0;
      mj=0;
      Vj=0;
      Gj=0;
      Qj=0;
      Ij=0;


      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Text(
              extent={{-146,-44},{154,-84}},
              textString="%name",
              lineColor={128,0,255})}),
        Documentation(revisions="<html>
<p><i>2009-2025</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ElectronTransfer;

    model ExternalSubstance "Constant source of molar concentration"
       extends Chemical.Undirected.Boundaries.Internal.PartialSubstanceInSolution(useSolution=false,
        useFore=false,
        useRear=false);

      parameter stateOfMatter.SubstanceDataParameters substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true, Dialog(enable=not useRear));

      parameter Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities quantity "Concentration quantity";

      parameter Real FixedValue = 1e-8
      "Fixed value of concentration in selected quantity if useVariableInput=false"
        annotation (HideResult=true, Dialog(enable=not useVariableInput));

      parameter Boolean useVariableInput = false
      "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      Modelica.Blocks.Interfaces.RealInput VariableInput=val if useVariableInput
        annotation (HideResult=true, Placement(transformation(extent={{-130,56},{-90,96}})));

      Real  value(unit=Chemical.Undirected.Boundaries.Internal.getUnit(quantity));
    equation
      if not useRear then
       substanceDataVar = substanceData;
      end if;

      if not useVariableInput then
        value=FixedValue;
      end if;

      //mole fraction
      //substance.x = val / solutionState.p;

      dH = 0;
      i = 0;
      dV = 0;
      Gj = 0;
      nj = 0;
      mj = 0;
      Qj = 0;
      Ij = 0;
      Vj = 0;

      if quantity == Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.c_molpm3 then
        value =  (substance.x * solutionState.n)/solutionState.V;
      elseif quantity == Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.X_kgpkg then
        value = ((substance.x * solutionState.n)/solutionState.m)/stateOfMatter.specificAmountOfParticles(substanceDataVar,
       solutionState.T,
       solutionState.p,
       solutionState.v,
       solutionState.I);
      elseif quantity == Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.b_molpkg then
        value =  (substance.x * solutionState.n)/solutionState.m;
      elseif quantity == Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.x_molpmol then
        value = substance.x;
      elseif quantity == Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_Pa then
        value =  substance.x*solutionState.p;
      elseif quantity == Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_kPa then
        value*1000 =  substance.x*solutionState.p;
      elseif quantity == Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg then
        value =  substance.x*solutionState.p * (760/101325);
      elseif quantity == Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_bar then
        value =  Modelica.Units.Conversions.to_bar(substance.x*solutionState.p);
      end if;

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
              lineColor={128,0,255}),
            Text(
              extent={{-104,-76},{100,-100}},
              lineColor={0,0,0},
              textString="%T K")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ExternalSubstance;

    model ExternalGas "Gas substance with defined partial pressure"
      extends Chemical.Undirected.Boundaries.Internal.PartialSubstanceInSolution(
        useFore=false,
        useRear=false,                            redeclare package stateOfMatter = gasModel);

       replaceable package gasModel = Chemical.Interfaces.IdealGasMSL constrainedby
        Chemical.Interfaces.StateOfMatter "Gas substance model"
        annotation (choices(
          choice(redeclare package gasModel =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package gasModel =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package gasModel =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

       parameter gasModel.SubstanceDataParameters substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true, Dialog(enable=not useRear));

      parameter Boolean usePartialPressureInput = false
      "=true, if fixed partial pressure is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.Units.SI.Pressure PartialPressure=1e-05
        "Fixed partial pressure if usePartialPressureInput=false" annotation (
         HideResult=true, Dialog(enable=not usePartialPressureInput));

      Modelica.Blocks.Interfaces.RealInput partialPressure(start=
            PartialPressure, final unit="Pa")=p if usePartialPressureInput
      "Partial pressure of gas = total pressure * gas fraction"
        annotation (HideResult=true,Placement(transformation(extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-100,72}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-100,72})));

      Modelica.Units.SI.Pressure p "Current partial pressure";

    equation
      if not useRear then
       substanceDataVar = substanceData;
      end if;

      if not usePartialPressureInput then
        p=PartialPressure;
      end if;

      //mole fraction
      substance.x = p / solutionState.p;

      dH = 0;
      i = 0;
      dV = 0;
      Gj = 0;
      nj = 0;
      mj = 0;
      Qj = 0;
      Ij = 0;
      Vj = 0;

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
              fillColor={159,159,223},
              fillPattern=FillPattern.Backward,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Text(
              extent={{54,108},{-46,8}},
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
              lineColor={128,0,255}),
            Text(
              extent={{-100,-102},{104,-126}},
              lineColor={0,0,0},
              textString="%T K")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ExternalGas;

    model BoundaryRear "Generic Boundary model (may act as source or sink)"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Chemical.Interfaces.SolutionStateParameters solutionState
        annotation (Dialog(enable=not solutionFromInput));
      parameter stateOfMatter.SubstanceDataParameters substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true);


      parameter Boolean solutionFromInput = false "Use input connector for solution?";
      parameter Boolean potentialFromInput = false "Use input connector for chemical potential";
      parameter Boolean enthalpyFromInput = false "Use input connector for molar enthalpy";

      parameter Modelica.Units.SI.MolarEnthalpy h0_par=0 "molar enthalpy set value"
        annotation (Dialog(enable=not enthalpyFromInput));
      parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "ChemicalPotential set value" annotation (Dialog(enable=not potentialFromInput));

      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
        annotation (Dialog(tab="Advanced"));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the boundary" annotation (Dialog(tab="Advanced"));

      Modelica.Blocks.Interfaces.RealInput u0_var(unit="J/mol")
                                                             if potentialFromInput "Chemical potential input connector [J/mol]"
        annotation (Placement(transformation(extent={{-40,40},{0,80}}), iconTransformation(extent={{-40,40},{0,80}})));
      Modelica.Blocks.Interfaces.RealInput h0_var(unit="J/mol")  if  enthalpyFromInput "Enthalpy input connector"
        annotation (Placement(transformation(extent={{-40,-40},{0,0}}), iconTransformation(extent={{-40,-20},{0,20}})));
      Interfaces.Fore fore(redeclare package stateOfMatter = stateOfMatter)
                           annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      Modelica.Units.SI.ChemicalPotential u_rearwards=fore.state_rearwards.u;

      Modelica.Blocks.Interfaces.RealInput u0(unit="J/mol") "Internal potential connector";
      Modelica.Blocks.Interfaces.RealInput h0(unit = "J/mol") "Internal enthalpy connector";

      Modelica.Units.SI.ChemicalPotential r;
      Chemical.Interfaces.SolutionState s "State of chemical solution";

    public
      Chemical.Interfaces.SolutionPort solution(T=s.T,p=s.p,v=s.v,n=s.n,m=s.m,V=s.V,G=s.G,Q=s.Q,I=s.I, i=0, dH=0, dV=0, nj=0, mj=0, Vj=0, Gj=0, Qj=0, Ij=0) if solutionFromInput
        annotation (Placement(transformation(extent={{-30,-70},{-10,-50}}), iconTransformation(extent={{-30,-70},{-10,-50}})));
    equation

      connect(u0_var, u0);
      if not potentialFromInput then
        u0 = u0_par;
      end if;

      connect(h0_var, h0);
      if not enthalpyFromInput then
         h0 = h0_par;
      end if;


      if not solutionFromInput then
        s=solutionState;
      end if;

      der(fore.n_flow)*L = fore.r-r;

      //if port.n_flow > 0 -> it is sink (r=u_set-u_in) else it is source (r=0)
      r = .Chemical.Undirected.Internal.regStep(fore.n_flow, u0 - u_rearwards, 0, n_flow_reg);

      fore.state_forwards = Chemical.Interfaces.SubstanceState(u=u0,h=h0);
      fore.solution=s;
      fore.definition=substanceData;


      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{0,76},{64,-84}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{0,80},{60,-80}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{60,0},{84,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{60,80},{60,-80}},
              color={158,66,200},
              thickness=0.5),
            Line(points={{44,80},{44,-80}}, color={255,255,255}),
            Line(
              points={{28,80},{28,-80}},
              color={255,255,255},
              thickness=0.5),
            Line(
              points={{12,80},{12,-80}},
              color={255,255,255},
              thickness=1)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>A undirected boundary that can act as source and sink, depending on the rest of the system. The Boundary_rear has to be connected to the rear end of your model and therefore has a fore port.</u>
<u>At positive massflow the fore port acts as an outlet and therefore the boundary_rear is a source.</u>
</html>"));
    end BoundaryRear;

    model TerminalRear "Fore Boundary that imposes n_flow = 0"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Chemical.Interfaces.SolutionStateParameters solutionState
        annotation (Dialog(enable=not solutionFromInput));
      parameter stateOfMatter.SubstanceDataParameters substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true);


      parameter Boolean solutionFromInput = false "Use input connector for solution?";

      parameter Modelica.Units.SI.Time TC=0.1 "Time constant for potential adaption" annotation (Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.MolarEnthalpy h=0 "Source enthalpy";

      parameter Modelica.Units.SI.ChemicalPotential u_0=0 "Initial potential";

      Chemical.Undirected.Interfaces.Fore fore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

    protected
      Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);
      Chemical.Interfaces.SolutionState s "State of chemical solution";

    public
      Chemical.Interfaces.SolutionPort solution(T=s.T,p=s.p,v=s.v,n=s.n,m=s.m,V=s.V,G=s.G,Q=s.Q,I=s.I, i=0, dH=0, dV=0, nj=0, mj=0, Vj=0, Gj=0, Qj=0, Ij=0) if solutionFromInput
        annotation (Placement(transformation(extent={{0,-30},{20,-10}}),    iconTransformation(extent={{0,-30},{20,-10}})));

    initial equation
      u = u_0;

    equation
      if not solutionFromInput then
        s=solutionState;
      end if;


      fore.n_flow = 0;

      TC * der(u) =fore.r;
      fore.state_forwards = Chemical.Interfaces.SubstanceState(u=u,h= h);


      fore.solution=s;
      fore.definition=substanceData;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{34,26},{74,-34}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{68,0},{84,0}},
              color={158,66,200},
              thickness=0.5),
            Rectangle(
              extent={{30,30},{70,-30}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{70,30},{30,-30}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{30,30},{70,-30}},
              color={158,66,200},
              thickness=0.5)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Rear Boundary that terminates the flow.  The Boundary has to be connected to the rear end of your model and therefore has a fore port.</u>
<u>It imposes a n_flow=0 boundary and with a time constant, adapts the potential sucht, that inertal potential r goes to zero.</u>
</html>"));
    end TerminalRear;

    model BoundaryFore "Generic Boundary model (may act as source or sink)"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));


      parameter Boolean potentialFromInput = false "Use input connector for potential?";
      parameter Boolean enthalpyFromInput = false "Use input connector for molar enthalpy";

      parameter Modelica.Units.SI.MolarEnthalpy h0_par=0 "molar enthalpy set value"
        annotation (Dialog(enable=not enthalpyFromInput));
      parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "ChemicalPotential set value" annotation (Dialog(enable=not potentialFromInput));

      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the boundary" annotation (Dialog(tab="Advanced"));

      Modelica.Blocks.Interfaces.RealInput u0_var(unit="J/mol") if potentialFromInput "ChemicalPotential input connector [Pa]"
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,60}),
          iconTransformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,60})));
      Modelica.Blocks.Interfaces.RealInput h0_var(unit = "J/mol") if enthalpyFromInput "Enthalpy input connector"
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,-20}),
          iconTransformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,0})));
      Chemical.Undirected.Interfaces.Rear rear(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-80,-20},{-120,20}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      Modelica.Units.SI.ChemicalPotential u_forwards=rear.state_forwards.u;

      Modelica.Blocks.Interfaces.RealInput u0(unit="J/mol") "Internal potential connector";
      Modelica.Blocks.Interfaces.RealInput h0(unit = "J/mol") "Internal enthalpy connector";

      Modelica.Units.SI.ChemicalPotential r;


    equation

      connect(u0_var, u0);
      if not potentialFromInput then
        u0 = u0_par;
      end if;

      connect(h0_var, h0);
      if not enthalpyFromInput then
         h0 = h0_par;
      end if;


      der(rear.n_flow)*L = rear.r-r;

      //if port.n_flow > 0 -> it is sink (r=u_set-u_in) else it is source (r=0)
      r = .Chemical.Undirected.Internal.regStep(rear.n_flow, u0 - u_forwards, 0, n_flow_reg);

      rear.state_rearwards = Chemical.Interfaces.SubstanceState(u=u0,h=h0);


      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{4,76},{-60,-84}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{0,80},{-60,-80}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{-60,0},{-84,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{-60,80},{-60,-80}},
              color={158,66,200},
              thickness=0.5),
            Line(points={{-44,80},{-44,-80}}, color={255,255,255}),
            Line(
              points={{-26,80},{-26,-80}},
              color={255,255,255},
              thickness=0.5),
            Line(
              points={{-12,80},{-12,-80}},
              color={255,255,255},
              thickness=1)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>A undirected boundary that can act as source and sink, depending on the rest of the system. The Boundary_fore has to be connected to the fore end of your model and therefore has a rear port.</u>
<u>At positive massflow the rear port acts as an inlet and therefore the boundary_fore is a sink.</u>
</html>"));
    end BoundaryFore;

    model TerminalFore "Rear Boundary that imposes n_flow = 0"
      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Modelica.Units.SI.Time TC=0.1 "Time constant for potential adaption" annotation (Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.MolarEnthalpy h=0 "Source enthalpy";

      parameter Modelica.Units.SI.ChemicalPotential u_0=0 "Initial potential";

      Chemical.Undirected.Interfaces.Rear rear(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

    protected
      Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

    initial equation
      u = u_0;

    equation
      rear.n_flow = 0;

      TC * der(u) =rear.r;
      rear.state_rearwards = Chemical.Interfaces.SubstanceState(u=u,h= h);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-20,30},{20,-30}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              origin={-46,-4},
              rotation=360),
            Line(
              points={{-10,0},{6,0}},
              color={158,66,200},
              thickness=0.5,
              origin={-74,0},
              rotation=360),
            Rectangle(
              extent={{-20,30},{20,-30}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              origin={-50,0},
              rotation=360),
            Line(
              points={{20,30},{-20,-30}},
              color={158,66,200},
              thickness=0.5,
              origin={-50,0},
              rotation=360),
            Line(
              points={{-20,30},{20,-30}},
              color={158,66,200},
              thickness=0.5,
              origin={-50,0},
              rotation=360)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Fore Boundary that terminates the flow. The boundary has to be connected to the fore end of your model and therefore has a rear port.</u>
<u>It imposes a n_flow=0 boundary and with a time constant, adapts the potential such, that inertial potential r goes to zero.</u>
</html>"));
    end TerminalFore;

    package Tests "Tests for the boundaries package"
      extends Modelica.Icons.ExamplesPackage;

      model TestSubstance
         extends Modelica.Icons.Example;
        Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,6}})));

        BoundaryRear boundaryRear(
          substanceData=Chemical.Substances.Water_liquid(),
          solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-76,14},{-56,34}})));
        Substance substance(
          useRear=true,
          useFore=false,
          substanceData=Chemical.Substances.Water_liquid())
                            annotation (Placement(transformation(extent={{24,14},{44,34}})));
        Substance substance2(
          useRear=false,
          useFore=true,
          useSolution=true,
          substanceData=Chemical.Substances.Water_liquid())
                                                  annotation (Placement(transformation(extent={{-68,-26},{-48,-6}})));
        BoundaryFore boundaryFore annotation (Placement(transformation(extent={{36,-26},{56,-6}})));
        Substance substance1(useRear=false,
          useFore=true,                     substanceData=Chemical.Substances.Water_liquid())
                                                  annotation (Placement(transformation(extent={{-72,72},{-52,92}})));
        BoundaryFore boundaryFore1
                                  annotation (Placement(transformation(extent={{32,72},{52,92}})));
        BoundaryRear boundaryRear1(substanceData=Chemical.Substances.Water_liquid(), solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-78,42},{-58,62}})));
        Substance substance3(
          useRear=true,
          useFore=true,      substanceData=Chemical.Substances.Water_liquid())
                            annotation (Placement(transformation(extent={{-28,42},{-8,62}})));
        BoundaryFore boundaryFore2
                                  annotation (Placement(transformation(extent={{30,42},{50,62}})));
        BoundaryRear boundaryRear2(substanceData=Chemical.Substances.Water_liquid(), solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-66,-82},{-46,-62}})));
        Substance substance4(
          useRear=true,
          useFore=false,
          useSolution=true,
          substanceData=Chemical.Substances.Water_liquid())
                            annotation (Placement(transformation(extent={{34,-82},{54,-62}})));
        BoundaryRear boundaryRear3(substanceData=Chemical.Substances.Water_liquid(), solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-68,-54},{-48,-34}})));
        Substance substance5(
          useRear=true,
          useFore=true,
          useSolution=true,                          substanceData=Chemical.Substances.Water_liquid())
                            annotation (Placement(transformation(extent={{-18,-54},{2,-34}})));
        BoundaryFore boundaryFore3
                                  annotation (Placement(transformation(extent={{40,-54},{60,-34}})));
      equation
        connect(boundaryRear.fore, substance.rear) annotation (Line(
            points={{-56,24},{24,24}},
            color={158,66,200},
            thickness=0.5));
        connect(substance2.fore,boundaryFore. rear) annotation (Line(
            points={{-48,-16},{36,-16}},
            color={158,66,200},
            thickness=0.5));
        connect(substance2.solution, solution.solution) annotation (Line(points={{-64,-26},{-64,-104},{60,-104},{60,-98.94}},        color={127,127,0}));
        connect(substance1.fore, boundaryFore1.rear) annotation (Line(
            points={{-52,82},{32,82}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear1.fore, substance3.rear) annotation (Line(
            points={{-58,52},{-28,52}},
            color={158,66,200},
            thickness=0.5));
        connect(substance3.fore, boundaryFore2.rear) annotation (Line(
            points={{-8,52},{30,52}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear2.fore, substance4.rear) annotation (Line(
            points={{-46,-72},{34,-72}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear3.fore, substance5.rear) annotation (Line(
            points={{-48,-44},{-18,-44}},
            color={158,66,200},
            thickness=0.5));
        connect(substance5.fore, boundaryFore3.rear) annotation (Line(
            points={{2,-44},{40,-44}},
            color={158,66,200},
            thickness=0.5));
        connect(substance5.solution, solution.solution) annotation (Line(points={{-14,-54},{-14,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
        connect(substance4.solution, solution.solution) annotation (Line(points={{38,-82},{38,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
         annotation (
          Icon(graphics,
               coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
          Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
      end TestSubstance;

      model TestElectronTransfer
         extends Modelica.Icons.Example;
        Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,0}})));
        Chemical.Undirected.Boundaries.ElectronTransfer electronTransfer(useRear=false,
          useFore=true,
          useSolution=true)                                                             annotation (Placement(transformation(extent={{-58,-32},{-38,-12}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore1(u0_par=1000)        annotation (Placement(transformation(extent={{56,-34},{76,-14}})));
        ElectronTransfer                                electronTransfer1(useRear=false, useFore=
              true)                                                                     annotation (Placement(transformation(extent={{-56,70},{-36,90}})));
        BoundaryFore                                boundaryFore2(u0_par=1000)        annotation (Placement(transformation(extent={{58,70},{78,90}})));
        ElectronTransfer electronTransfer2(useRear=true, useFore=true)
                                           annotation (Placement(transformation(extent={{-4,46},{16,66}})));
        ElectronTransfer electronTransfer3(useRear=true,
                                           useFore=false) annotation (Placement(transformation(extent={{50,20},{70,40}})));
        BoundaryRear boundaryRear(substanceData=Chemical.Substances.Electrone_solid()) annotation (Placement(transformation(extent={{-60,46},{-40,66}})));
        BoundaryRear boundaryRear1(substanceData=Chemical.Substances.Electrone_solid()) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        BoundaryFore boundaryFore annotation (Placement(transformation(extent={{58,46},{78,66}})));
        ElectronTransfer electronTransfer4(
          useRear=true,
          useFore=true,
          useSolution=true)                                        annotation (Placement(transformation(extent={{0,-62},{20,-42}})));
        ElectronTransfer electronTransfer5(
          useRear=true,                    useFore=false,
          useSolution=true)                                                       annotation (Placement(transformation(extent={{54,-88},{74,-68}})));
        BoundaryRear boundaryRear2(substanceData=Chemical.Substances.Electrone_solid()) annotation (Placement(transformation(extent={{-56,-62},{-36,-42}})));
        BoundaryRear boundaryRear3(substanceData=Chemical.Substances.Electrone_solid()) annotation (Placement(transformation(extent={{-56,-88},{-36,-68}})));
        BoundaryFore boundaryFore3 annotation (Placement(transformation(extent={{62,-62},{82,-42}})));
      equation
        connect(electronTransfer.fore, boundaryFore1.rear) annotation (Line(
            points={{-38,-22},{56,-22},{56,-24}},
            color={158,66,200},
            thickness=0.5));
        connect(electronTransfer.solution, solution.solution) annotation (Line(points={{-54,-32},{-54,-106},{60,-106},{60,-99}}, color={127,127,0}));
        connect(electronTransfer1.fore, boundaryFore2.rear) annotation (Line(
            points={{-36,80},{58,80}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear.fore, electronTransfer2.rear) annotation (Line(
            points={{-40,56},{-4,56}},
            color={158,66,200},
            thickness=0.5));
        connect(electronTransfer2.fore, boundaryFore.rear) annotation (Line(
            points={{16,56},{58,56}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear1.fore, electronTransfer3.rear) annotation (Line(
            points={{-40,30},{50,30}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear2.fore, electronTransfer4.rear) annotation (Line(
            points={{-36,-52},{0,-52}},
            color={158,66,200},
            thickness=0.5));
        connect(electronTransfer4.fore, boundaryFore3.rear) annotation (Line(
            points={{20,-52},{62,-52}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear3.fore, electronTransfer5.rear) annotation (Line(
            points={{-36,-78},{54,-78}},
            color={158,66,200},
            thickness=0.5));
        connect(electronTransfer4.solution, solution.solution) annotation (Line(points={{4,-62},{6,-62},{6,-106},{60,-106},{60,-99}}, color={127,127,0}));
        connect(electronTransfer5.solution, solution.solution) annotation (Line(points={{58,-88},{58,-94},{60,-94},{60,-99}}, color={127,127,0}));
         annotation (
          Icon(graphics,
               coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
          Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
      end TestElectronTransfer;

      model TestExternalSubstance
         extends Modelica.Icons.Example;
        Chemical.Solution solution(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasMSL                                                        "Ideal Gas from MSL")
                                   annotation (Placement(transformation(extent={{-100,-100},{100,6}})));

        replaceable package gasModel = Chemical.Interfaces.IdealGasMSL constrainedby
          Chemical.Interfaces.StateOfMatter "Gas substance model"
          annotation (choices(
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGas        "Ideal Gas"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

        ExternalSubstance
                  externalSubstance1(
          useRear=true,
          useFore=false,
          redeclare package stateOfMatter = gasModel,
          quantity=Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
          FixedValue=1)     annotation (Placement(transformation(extent={{24,14},{44,34}})));

        ExternalSubstance
                  externalSubstance2(
          useRear=false,
          useFore=true,
          useSolution=true,
          redeclare package stateOfMatter = gasModel,
          substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
          quantity=Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
          FixedValue=10)                          annotation (Placement(transformation(extent={{-68,-26},{-48,-6}})));

        BoundaryFore boundaryFore(redeclare package stateOfMatter = gasModel)
                                  annotation (Placement(transformation(extent={{36,-26},{56,-6}})));
        ExternalSubstance externalIdealGas(
          useRear=false,
          useFore=true,
          quantity=Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
          redeclare package stateOfMatter = gasModel,
          substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
          FixedValue=15)                                        annotation (Placement(transformation(extent={{-72,72},{-52,92}})));
        BoundaryFore boundaryFore1(redeclare package stateOfMatter = gasModel)
                                  annotation (Placement(transformation(extent={{32,72},{52,92}})));
        BoundaryRear boundaryRear1(
          redeclare package stateOfMatter = gasModel,
          substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
          solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-76,42},{-56,62}})));
        ExternalSubstance
                  externalSubstance(
          useRear=true,
          useFore=true,       redeclare package stateOfMatter = gasModel,
          quantity=Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
          FixedValue=10)    annotation (Placement(transformation(extent={{-26,42},{-6,62}})));

        BoundaryFore boundaryFore2(redeclare package stateOfMatter = gasModel)
                                  annotation (Placement(transformation(extent={{30,42},{50,62}})));
        BoundaryRear boundaryRear2(
          redeclare package stateOfMatter = gasModel,
                                   substanceData=Chemical.Substances.IdealGasesMSL.H2O(), solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-66,-82},{-46,-62}})));
        ExternalSubstance
                  externalSubstance4(
          useRear=true,        useFore=false,
          useSolution=true,
          redeclare package stateOfMatter = gasModel,
          quantity=Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
          FixedValue=1)     annotation (Placement(transformation(extent={{34,-82},{54,-62}})));

        BoundaryRear boundaryRear3(
          redeclare package stateOfMatter = gasModel,
                                   substanceData=Chemical.Substances.IdealGasesMSL.H2O(), solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-68,-54},{-48,-34}})));
        ExternalSubstance
                  externalSubstance3(
          useRear=true,
          useFore=true,
          useSolution=true,
          redeclare package stateOfMatter = gasModel,
          quantity=Chemical.Undirected.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
          FixedValue=1)     annotation (Placement(transformation(extent={{-18,-54},{2,-34}})));

        BoundaryFore boundaryFore3(redeclare package stateOfMatter = gasModel)
                                  annotation (Placement(transformation(extent={{40,-54},{60,-34}})));
        BoundaryRear boundaryRear4(
          redeclare package stateOfMatter = gasModel,
          substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
          solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-74,14},{-54,34}})));
      equation
        connect(externalSubstance2.fore, boundaryFore.rear) annotation (Line(
            points={{-48,-16},{36,-16}},
            color={158,66,200},
            thickness=0.5));
        connect(externalSubstance2.solution, solution.solution) annotation (Line(points={{-64,-26},{-64,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
        connect(externalIdealGas.fore, boundaryFore1.rear) annotation (Line(
            points={{-52,82},{32,82}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear1.fore, externalSubstance.rear) annotation (Line(
            points={{-56,52},{-26,52}},
            color={158,66,200},
            thickness=0.5));
        connect(externalSubstance.fore, boundaryFore2.rear) annotation (Line(
            points={{-6,52},{30,52}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear2.fore, externalSubstance4.rear) annotation (Line(
            points={{-46,-72},{34,-72}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear3.fore, externalSubstance3.rear) annotation (Line(
            points={{-48,-44},{-18,-44}},
            color={158,66,200},
            thickness=0.5));
        connect(externalSubstance3.fore, boundaryFore3.rear) annotation (Line(
            points={{2,-44},{40,-44}},
            color={158,66,200},
            thickness=0.5));
        connect(externalSubstance3.solution, solution.solution) annotation (Line(points={{-14,-54},{-14,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
        connect(externalSubstance4.solution, solution.solution) annotation (Line(points={{38,-82},{38,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
        connect(boundaryRear4.fore, externalSubstance1.rear) annotation (Line(
            points={{-54,24},{24,24}},
            color={158,66,200},
            thickness=0.5));
         annotation (
          Icon(graphics,
               coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
          Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
      end TestExternalSubstance;

      model TestExternalGas
         extends Modelica.Icons.Example;
        Chemical.Solution solution(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasMSL                                                        "Ideal Gas from MSL")
                                   annotation (Placement(transformation(extent={{-100,-100},{100,6}})));

        replaceable package gasModel = Chemical.Interfaces.IdealGasMSL constrainedby
          Chemical.Interfaces.StateOfMatter "Gas substance model"
          annotation (choices(
            choice(redeclare package gasModel =
              Chemical.Interfaces.IdealGas        "Ideal Gas"),
            choice(redeclare package gasModel =
              Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
            choice(redeclare package gasModel =
              Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

        ExternalGas
                  externalGas1(
          useRear=true,
          useFore=false,
          redeclare package gasModel = gasModel,
                         PartialPressure(displayUnit="mmHg") = 133.322387415)
                            annotation (Placement(transformation(extent={{24,14},{44,34}})));
        ExternalGas
                  externalGas2(
          useRear=false,
          useFore=true,
          useSolution=true,
          redeclare package gasModel = gasModel,
          substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
          PartialPressure(displayUnit="mmHg") = 1333.22387415)
                                                  annotation (Placement(transformation(extent={{-68,-26},{-48,-6}})));
        BoundaryFore boundaryFore(redeclare package stateOfMatter = gasModel)
                                  annotation (Placement(transformation(extent={{36,-26},{56,-6}})));
        ExternalGas externalIdealGas(
          useRear=false,
          useFore=true,
          redeclare package gasModel = gasModel,
          substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
          PartialPressure(displayUnit="mmHg") = 1999.835811225) annotation (Placement(transformation(extent={{-72,72},{-52,92}})));
        BoundaryFore boundaryFore1(redeclare package stateOfMatter = gasModel)
                                  annotation (Placement(transformation(extent={{32,72},{52,92}})));
        BoundaryRear boundaryRear1(
          redeclare package stateOfMatter = gasModel,
          substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
          solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-76,42},{-56,62}})));
        ExternalGas
                  externalGas(
          useRear=true,
          useFore=true,       redeclare package gasModel = gasModel,                 PartialPressure(displayUnit="mmHg") = 1333.22387415)
                            annotation (Placement(transformation(extent={{-26,42},{-6,62}})));
        BoundaryFore boundaryFore2(redeclare package stateOfMatter = gasModel)
                                  annotation (Placement(transformation(extent={{30,42},{50,62}})));
        BoundaryRear boundaryRear2(
          redeclare package stateOfMatter = gasModel,
                                   substanceData=Chemical.Substances.IdealGasesMSL.H2O(), solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-66,-82},{-46,-62}})));
        ExternalGas
                  externalGas4(
          useRear=true,        useFore=false,
          useSolution=true,
          redeclare package gasModel = gasModel,
          PartialPressure(displayUnit="mmHg") = 133.322387415)
                            annotation (Placement(transformation(extent={{34,-82},{54,-62}})));
        BoundaryRear boundaryRear3(
          redeclare package stateOfMatter = gasModel,
                                   substanceData=Chemical.Substances.IdealGasesMSL.H2O(), solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-68,-54},{-48,-34}})));
        ExternalGas
                  externalGas3(
          useRear=true,
          useFore=true,
          useSolution=true,
          redeclare package gasModel = gasModel,       PartialPressure(displayUnit="mmHg") = 133.322387415)
                            annotation (Placement(transformation(extent={{-18,-54},{2,-34}})));
        BoundaryFore boundaryFore3(redeclare package stateOfMatter = gasModel)
                                  annotation (Placement(transformation(extent={{40,-54},{60,-34}})));
        BoundaryRear boundaryRear4(
          redeclare package stateOfMatter = gasModel,
          substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
          solutionFromInput=false)
                         annotation (Placement(transformation(extent={{-74,14},{-54,34}})));
      equation
        connect(externalGas2.fore, boundaryFore.rear) annotation (Line(
            points={{-48,-16},{36,-16}},
            color={158,66,200},
            thickness=0.5));
        connect(externalGas2.solution, solution.solution) annotation (Line(points={{-64,-26},{-64,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
        connect(externalIdealGas.fore, boundaryFore1.rear) annotation (Line(
            points={{-52,82},{32,82}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear1.fore, externalGas.rear) annotation (Line(
            points={{-56,52},{-26,52}},
            color={158,66,200},
            thickness=0.5));
        connect(externalGas.fore, boundaryFore2.rear) annotation (Line(
            points={{-6,52},{30,52}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear2.fore, externalGas4.rear) annotation (Line(
            points={{-46,-72},{34,-72}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear3.fore, externalGas3.rear) annotation (Line(
            points={{-48,-44},{-18,-44}},
            color={158,66,200},
            thickness=0.5));
        connect(externalGas3.fore, boundaryFore3.rear) annotation (Line(
            points={{2,-44},{40,-44}},
            color={158,66,200},
            thickness=0.5));
        connect(externalGas3.solution, solution.solution) annotation (Line(points={{-14,-54},{-14,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
        connect(externalGas4.solution, solution.solution) annotation (Line(points={{38,-82},{38,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
        connect(boundaryRear4.fore, externalGas1.rear) annotation (Line(
            points={{-54,24},{24,24}},
            color={158,66,200},
            thickness=0.5));
         annotation (
          Icon(graphics,
               coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
          Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
      end TestExternalGas;

      model TestBoundaries "Tests for the rear and fore boundary"
        extends Modelica.Icons.Example;

        BoundaryRear boundary_rear(
          u0_par=100000,
          fore(n_flow(start=0, fixed=true)))
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-28,82})));
        BoundaryFore boundary_fore(
          potentialFromInput=true,
          u0_par=110000) annotation (Placement(transformation(extent={{22,72},{42,92}})));
        inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=0.01) annotation (Placement(transformation(extent={{-88,72},{-68,92}})));
        Modelica.Blocks.Sources.Step step(
          height=-100000,
          offset=140000,
          startTime=5)
          annotation (Placement(transformation(extent={{62,76},{50,88}})));
        TerminalRear terminal_rear(h=0, u_0=0)
          annotation (Placement(transformation(extent={{-38,46},{-18,66}})));
        BoundaryFore boundary_fore1(
          potentialFromInput=true,
          u0_par=110000) annotation (Placement(transformation(extent={{22,46},{42,66}})));
        Modelica.Blocks.Sources.Step step1(
          height=-100000,
          offset=140000,
          startTime=5)
          annotation (Placement(transformation(extent={{62,50},{50,62}})));
        BoundaryRear boundary_rear1(
          potentialFromInput=true,
          u0_par=100000) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-26,30})));
        TerminalFore terminal_fore(h=0, u_0=0)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={34,30})));
        Modelica.Blocks.Sources.Step step2(
          height=-100000,
          offset=140000,
          startTime=5)
          annotation (Placement(transformation(extent={{-52,24},{-40,36}})));
        BoundaryRear boundary_rear2(
          solutionFromInput=true,
          u0_par=100000,
          fore(n_flow(start=0, fixed=true)))
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-60,-22})));
        BoundaryFore boundary_fore2(potentialFromInput=true, u0_par=110000)
                         annotation (Placement(transformation(extent={{-10,-32},{10,-12}})));
        Modelica.Blocks.Sources.Step step3(
          height=-100000,
          offset=140000,
          startTime=5)
          annotation (Placement(transformation(extent={{30,-28},{18,-16}})));
        TerminalRear terminal_rear1(solutionFromInput=true)
          annotation (Placement(transformation(extent={{-36,-60},{-16,-40}})));
        BoundaryFore boundary_fore3(potentialFromInput=true, u0_par=110000)
                         annotation (Placement(transformation(extent={{24,-60},{44,-40}})));
        Modelica.Blocks.Sources.Step step4(
          height=-100000,
          offset=140000,
          startTime=5)
          annotation (Placement(transformation(extent={{64,-56},{52,-44}})));
        BoundaryRear boundary_rear3(
          substanceData=Chemical.Substances.Water_liquid(),
          solutionFromInput=true,
          potentialFromInput=true,
          u0_par=100000) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={26,-78})));
        TerminalFore terminal_fore1
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={82,-78})));
        Modelica.Blocks.Sources.Step step5(
          height=-100000,
          offset=140000,
          startTime=5)
          annotation (Placement(transformation(extent={{-10,-90},{2,-78}})));
        Solution solution annotation (Placement(transformation(extent={{-98,-98},{102,0}})));
      equation
        connect(step.y, boundary_fore.u0_var)
          annotation (Line(points={{49.4,82},{42,82},{42,88},{34,88}},
                                                         color={0,0,127}));
        connect(boundary_fore.rear, boundary_rear.fore) annotation (Line(
            points={{22,82},{-18,82}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore1.rear, terminal_rear.fore) annotation (Line(
            points={{22,56},{-18,56}},
            color={158,66,200},
            thickness=0.5));
        connect(step1.y, boundary_fore1.u0_var) annotation (Line(points={{49.4,56},{42,56},{42,62},{34,62}},
                                                                                               color={0,0,127}));
        connect(step2.y,boundary_rear1.u0_var)  annotation (Line(points={{-39.4,30},{-34,30},{-34,24},{-28,24}},
                                                                                                 color={0,0,127}));
        connect(boundary_rear1.fore, terminal_fore.rear) annotation (Line(
            points={{-16,30},{24,30}},
            color={158,66,200},
            thickness=0.5));
        connect(step3.y, boundary_fore2.u0_var) annotation (Line(points={{17.4,-22},{10,-22},{10,-16},{2,-16}}, color={0,0,127}));
        connect(boundary_fore2.rear, boundary_rear2.fore) annotation (Line(
            points={{-10,-22},{-50,-22}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore3.rear, terminal_rear1.fore) annotation (Line(
            points={{24,-50},{-16,-50}},
            color={158,66,200},
            thickness=0.5));
        connect(step4.y,boundary_fore3. u0_var) annotation (Line(points={{51.4,-50},{44,-50},{44,-44},{36,-44}},
                                                                                               color={0,0,127}));
        connect(step5.y,boundary_rear3.u0_var)  annotation (Line(points={{2.6,-84},{24,-84}},    color={0,0,127}));
        connect(boundary_rear3.fore, terminal_fore1.rear) annotation (Line(
            points={{36,-78},{72,-78}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear2.solution, solution.solution) annotation (Line(points={{-62,-16},{-62,-102},{62,-102},{62,-97.02}}, color={127,127,0}));
        connect(terminal_rear1.solution, solution.solution) annotation (Line(points={{-25,-52},{-26,-52},{-26,-102},{62,-102},{62,-97.02}}, color={127,127,0}));
        connect(boundary_rear3.solution,solution. solution) annotation (Line(points={{24,-72},{24,-102},{62,-102},{62,-97.02}},
                                                                                                                    color={127,127,0}));
        annotation (
          Icon(graphics,
               coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
          Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
      end TestBoundaries;
      annotation (Documentation(info="<html>
<u>Tests for the boundaries package.</u>
</html>"));
    end Tests;

    package Internal "Partials and Internal functions"
    extends Modelica.Icons.InternalPackage;

      partial model PartialSubstance

        replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
          Chemical.Interfaces.StateOfMatter
        "Substance model to translate data into substance properties"
          annotation (choices(
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.Incompressible  "Incompressible"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGas        "Ideal Gas"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

        outer Modelica.Fluid.System system "System wide properties";

        parameter Boolean useRear = true "Use rearward conector?"
            annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));
        parameter Boolean useFore = true "Use forward connector?"
            annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

        parameter Boolean initialize_potential = true "If true: initialize ChemicalPotential"
          annotation(Dialog(tab= "Initialization"));
        parameter Modelica.Units.SI.ChemicalPotential u_start=0 "Initial ChemicalPotential" annotation (Dialog(tab="Initialization", enable=initialize_potential));

        parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance at inlet and outlet" annotation (Dialog(tab="Advanced"));
        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
          annotation (Dialog(tab="Advanced"));

        stateOfMatter.BaseProperties substance;

        Chemical.Undirected.Interfaces.Rear rear(
          redeclare package stateOfMatter = stateOfMatter,
          n_flow=n_flow_rear,
          r=r_rear_port,
          state_rearwards=state_out_rear,
          state_forwards=state_in_rear,
          solution=solutionState,
          definition=substanceDataVar) if useRear
          annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

        Chemical.Undirected.Interfaces.Fore fore(
          redeclare package stateOfMatter = stateOfMatter,
          n_flow=n_flow_fore,
          r=r_fore_port,
          state_forwards=state_out_fore,
          state_rearwards=state_in_fore,
          solution=solutionState,
          definition=substanceDataVar) if useFore
          annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

         //Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential of the substance";
         //Modelica.Units.SI.MolarEnthalpy molarEnthalpy "Molar enthalpy of the substance";
         Modelica.Units.SI.MolarFlowRate n_flow "Molar change of the amount of base substance";
         Modelica.Units.SI.EnthalpyFlowRate h_flow "Change of enthalpy";

       protected
        outer Chemical.DropOfCommons dropOfCommons;

        stateOfMatter.SubstanceData substanceDataVar;
        Chemical.Interfaces.SolutionState solutionState;

         //if port.n_flow > 0 -> it is sink (r=medium.u-u_in) else it is source (r=0)
        Modelica.Units.SI.ChemicalPotential r_rear_intern=Chemical.Undirected.Internal.regStep(
                  n_flow_rear,
                  substance.u - state_in_rear.u,
                  0,
                  n_flow_reg);
        Modelica.Units.SI.ChemicalPotential r_fore_intern=Chemical.Undirected.Internal.regStep(
                  n_flow_fore,
                  substance.u - state_in_fore.u,
                  0,
                  n_flow_reg);
        // dont regstep variables that are only in der(state), to increase accuracy
        Modelica.Units.SI.EnthalpyFlowRate h_flow_rear=(if n_flow_rear >= 0 then state_in_rear.h else h_out_rear)*n_flow_rear;
        Modelica.Units.SI.EnthalpyFlowRate h_flow_fore=(if n_flow_fore >= 0 then state_in_fore.h else h_out_fore)*n_flow_fore;

        Chemical.Interfaces.SubstanceState state_out_rear;
        Modelica.Units.SI.MolarEnthalpy h_out_rear=state_out_rear.h;

        Chemical.Interfaces.SubstanceState state_out_fore;
        Modelica.Units.SI.MolarEnthalpy h_out_fore=state_out_fore.h;

        Modelica.Units.SI.ChemicalPotential r_rear_port;
        Modelica.Units.SI.ChemicalPotential r_fore_port;
        Modelica.Units.SI.MolarFlowRate n_flow_rear;
        Modelica.Units.SI.MolarFlowRate n_flow_fore;

        InputSubstanceData state_in_rear;
        InputSubstanceData state_in_fore;

        connector InputSubstanceData = input Chemical.Interfaces.SubstanceState
          "Substance definition as input signal connector";
      equation
        substance.substanceDataVar=substanceDataVar;

        substance.T = solutionState.T;
        substance.p = solutionState.p;
        substance.v = solutionState.v;
        substance.I = solutionState.I;


        state_out_rear = Chemical.Interfaces.SubstanceState(u=substance.u,h=substance.h);
        state_out_fore = state_out_rear;


        der(n_flow_rear)*L = r_rear_port - r_rear_intern;
        der(n_flow_fore)*L = r_fore_port - r_fore_intern;

        n_flow = n_flow_rear + n_flow_fore;
        if not useRear then
          h_flow = h_flow_fore;
        elseif not useFore then
          h_flow = h_flow_rear;
        else
          h_flow = h_flow_rear + h_flow_fore;
        end if;


        if not useRear then
          n_flow_rear = 0;
          state_in_rear = Chemical.Interfaces.SubstanceState(u=substance.u,h=substance.h);
          end if;
        if not useFore then
          n_flow_fore = 0;
          state_in_fore = Chemical.Interfaces.SubstanceState(u=substance.u,h=substance.h);
        end if;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end PartialSubstance;

      partial model PartialSubstanceInSolution "Substance properties for components, where the substance is connected with the solution"
       extends PartialSubstance;

        parameter Boolean useSolution = false "Use solution connector?"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

        parameter Chemical.Interfaces.SolutionStateParameters solutionParam "Constant chemical solution state if not from rear or input"
          annotation (Dialog(enable=not useSolution and not useRear));

        Chemical.Interfaces.SolutionPort solution(
            T=solutionPortState.T,
            p=solutionPortState.p,
            v=solutionPortState.v,
            n=solutionPortState.n,
            m=solutionPortState.m,
            V=solutionPortState.V,
            G=solutionPortState.G,
            Q=solutionPortState.Q,
            I=solutionPortState.I,
            i=i,
            dH=dH,
            dV=dV,
            nj=nj,
            mj=mj,
            Vj=Vj,
            Gj=Gj,
            Qj=Qj,
            Ij=Ij)
              if useSolution "To connect substance with solution, where is pressented"
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

      protected
           Real i,dH,dV,nj,mj,Vj,Gj,Qj,Ij;

           Chemical.Interfaces.SolutionState solutionPortState;
      equation


        if (useSolution and not useRear) or (not useSolution) then
          solutionState=solutionPortState;
        end if;

        if not useSolution and not useRear then
          solutionState.T=solutionParam.T "Temperature of the solution";
          solutionState.p=solutionParam.p "Pressure of the solution";
          solutionState.v=solutionParam.v "Electric potential in the solution";
          solutionState.n=solutionParam.n "Amount of the solution";
          solutionState.m=solutionParam.m "Mass of the solution";
          solutionState.V=solutionParam.V "Volume of the solution";
          solutionState.G=solutionParam.G "Free Gibbs energy of the solution";
          solutionState.Q=solutionParam.Q "Electric charge of the solution";
          solutionState.I=solutionParam.I "Mole fraction based ionic strength of the solution";
        end if;

      end PartialSubstanceInSolution;

      partial model ConditionalSolutionFlow "Input of solution molar flow vs. parametric solution molar flow"

        parameter Boolean useSolutionFlowInput = false
        "=true, if solution flow is provided via input"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),
                Dialog(group="Conditional inputs", __Dymola_compact=true));

      parameter Modelica.Units.SI.VolumeFlowRate SolutionFlow=0
        "Volume flow rate of the solution if useSolutionFlowInput=false"
        annotation (HideResult=true, Dialog(enable=not useSolutionFlowInput));

        Modelica.Blocks.Interfaces.RealInput solutionFlow(start=SolutionFlow, final unit="m3/s")=
           volumeFlow if useSolutionFlowInput
           annotation ( HideResult=true, Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=270,
              origin={0,40}), iconTransformation(
              extent={{-20,-20},{20,20}},
              rotation=270,
              origin={0,40})));

      Modelica.Units.SI.VolumeFlowRate volumeFlow "Current solution volume flow";

      equation
        if not useSolutionFlowInput then
          volumeFlow = SolutionFlow;
        end if;

      end ConditionalSolutionFlow;

      package Types
        type ConcentrationQuantities = enumeration(
            c_molpm3 "Concentration (mmol/L)",
            X_kgpkg "Mass fraction (kg/kg)",
            b_molpkg "Molality (mol/kg)",
            x_molpmol "Mole fraction (mol/mol)",
            p_Pa "Partial pressure (Pa)",
            p_kPa "Partial pressure (kPa)",
            p_mmHg "Partial pressure (mmHg)",
            p_bar "Partial pressure (bar)");
      end Types;

      function getUnit "Returns unit of input quantity"
        extends Modelica.Icons.Function;

        input Types.ConcentrationQuantities quantity;
        output String unit;

      algorithm

        if quantity == Types.ConcentrationQuantities.c_molpm3 then
          unit := "mol/m3";
        elseif quantity == Types.ConcentrationQuantities.X_kgpkg then
          unit := "kg/kg";
        elseif quantity == Types.ConcentrationQuantities.b_molpkg then
          unit := "mol/kg";
        elseif quantity == Types.ConcentrationQuantities.x_molpmol then
          unit := "mol/mol";
        elseif quantity == Types.ConcentrationQuantities.p_Pa then
          unit := "Pa";
        elseif quantity == Types.ConcentrationQuantities.p_kPa then
          unit := "kPa";
        elseif quantity == Types.ConcentrationQuantities.p_mmHg then
          unit := "mmHg";
        elseif quantity == Types.ConcentrationQuantities.p_bar then
          unit := "bar";
        else
          unit :="";
        end if;

        annotation (Documentation(info="<html>
<p>Helper function to get the unit for a quantity.</p>
</html>"));
      end getUnit;
    annotation (Documentation(info="<html>
<p>This package contains all internal functions, partials, and other (e.g. experimental) models for the Boundaries package.</p>
</html>"));
    end Internal;
  annotation (Documentation(info="<html>
<u>
Boundary models for undirected chemical simulation.
</u>
</html>"));
  end Boundaries;

  package Topology "Junctions and Connectors for undirected chemical simulation"
    extends Modelica.Icons.Package;

    model JunctionMN "Generalized junction/splitter for undirected flow"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));


      parameter Integer N(min=0) = 1 "Number of rears";
      parameter Integer M(min=0) = 1 "Number of fors";

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance for each branch" annotation (Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of molar flow rate"
        annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Rear rears[N](redeclare package stateOfMatter =
            stateOfMatter)                                                                          "Rear ports"
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
      Chemical.Undirected.Interfaces.Fore fores[M](redeclare package stateOfMatter =
            stateOfMatter)                                                                          "Fore ports"
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={100,0})));

      Modelica.Units.SI.ChemicalPotential u_mix "mixing u assuming positive molarFlow";

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      Modelica.Units.SI.ChemicalPotential rs[M + N] "rs of the connectors internal of L";
      Modelica.Units.SI.ChemicalPotential r_mix "mixed r assuming positive molarflow";

      Modelica.Units.SI.MolarFlowRate inflows[M + N] "molarFlows assuming positive molarFlow";

      Modelica.Units.SI.ChemicalPotential us[M + N] "potential of medium entering";
      Modelica.Units.SI.ChemicalPotential us_out[M + N] "potential of medium exeting";

      Modelica.Units.SI.MolarEnthalpy hs[M + N] "molar enthalpy of medium entering";
      Modelica.Units.SI.MolarEnthalpy hs_out[M + N] "molar enthalpy of medium exeting";




    equation
      // rears are 1:N
      der(rears.n_flow)*L = rears.r-rs[1:N];
      for i in 1:N loop
         // inputs ports
        us[i] = rears[i].state_forwards.u;
        hs[i] = rears[i].state_forwards.h;
        // inflows, u_equation and set output
        inflows[i] = max(rears[i].n_flow, n_flow_reg);
        Chemical.Undirected.Internal.regStep(
              rears[i].n_flow,
              us[i],
              us_out[i],
              n_flow_reg) + rs[i] = u_mix + r_mix;
        rears[i].state_rearwards = Chemical.Interfaces.SubstanceState(u=us_out[i],h=hs_out[i]);
      end for;

      // fores are N+1:end
      der(fores.n_flow)*L = fores.r-rs[N+1:end];
      for i in 1:M loop
         // inputs ports
        us[N+i] = fores[i].state_rearwards.u;
        hs[N+i] = fores[i].state_rearwards.h;
        // inflows, u_equation and set output
        inflows[N+i] = max(fores[i].n_flow, n_flow_reg);
        Chemical.Undirected.Internal.regStep(
              fores[i].n_flow,
              us[N + i],
              us_out[N + i],
              n_flow_reg) + rs[N + i] = u_mix + r_mix;
        fores[i].state_forwards = Chemical.Interfaces.SubstanceState(u=us_out[N+i],h=hs_out[N+i]);

        fores[i].definition = rears[1].definition;
        fores[i].solution = rears[1].solution;

      end for;

      // Mass balance
      sum(rears.n_flow) + sum(fores.n_flow) = 0;
      //compute u_mix for rs computation
      u_mix =(us*inflows)/sum(inflows);

      // compute output quantities
      for i in 1:M+N loop
        hs_out[i] = (hs*inflows - hs[i]*inflows[i]) /(sum(inflows) - inflows[i]);
        us_out[i] = (us*inflows - us[i]*inflows[i]) /(sum(inflows) - inflows[i]);
      end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-84,0},{84,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{60,20},{100,60}},
              textColor={175,175,175},
              textString="%M"),
            Text(
              extent={{-100,20},{-60,60}},
              textColor={175,175,175},
              textString="%N")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>This model represents a generalized junction/splitter for undirected flow with N rear and M fore ports. </u>
<u>Note that in the undirected case a distinction between junction and splitter is not possible, since the flow direction is unknown in advance. </u>
</html>"));
    end JunctionMN;

    model ConnectInletFore
      "Directed/undirected connector with input and fore"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Inlet inlet(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}),
            iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
      Chemical.Undirected.Interfaces.Fore fore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      Modelica.Units.SI.ChemicalPotential r_fore;
      Modelica.Units.SI.ChemicalPotential r_inlet;

    equation
      fore.state_forwards = inlet.state;

      fore.n_flow + inlet.n_flow = 0;
      r_inlet + inlet.state.u = r_fore + Chemical.Undirected.Internal.regStep(
            inlet.n_flow,
            inlet.state.u,
            fore.state_rearwards.u,
            dropOfCommons.n_flow_reg);

      L/2*der(inlet.n_flow) = inlet.r - r_inlet;
      L/2*der(fore.n_flow) = fore.r - r_fore;

      fore.definition = inlet.definition;
      fore.solution = inlet.solution;

      annotation (Icon(
          graphics={
            Line(
              points={{-30,0},{30,0}},
              color={158,66,200},
              thickness=0.5)},
          coordinateSystem(preserveAspectRatio=false)),
         Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>This connector can be used to connect a unidirectional outlet to a undirected rear port. </u>
<u>The state from the inlet is given to the forward direction of the fore port, the total potential as well as the massflow of inlet and port are set equal. </u>
<u>Note that when the flow is reversed, the resulting inertial potential can be different on both sides of this connector. </u>
</html>"));
    end ConnectInletFore;

    model ConnectInletRear
      "Directed/undirected connector with input and rear"

     replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Inlet inlet(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-40,0}),
            iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
      Chemical.Undirected.Interfaces.Rear rear(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={40,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));
      ConnectInletFore connectInletFore(redeclare package stateOfMatter =
            stateOfMatter, final L=L/2)
        annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
      ConnectRearRear connectRearRear(redeclare package stateOfMatter =
            stateOfMatter, final L=L/2)
        annotation (Placement(transformation(extent={{2,-10},{22,10}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(connectInletFore.fore, connectRearRear.rear_a) annotation (Line(
          points={{-7,0},{9,0}},
          color={158,66,200},
          thickness=0.5));
      connect(inlet, connectInletFore.inlet) annotation (Line(
          points={{-40,0},{-13,0}},
          color={158,66,200},
          thickness=0.5));
      connect(connectRearRear.rear_b, rear) annotation (Line(
          points={{15,0},{40,0}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(
          graphics={
            Line(
              points={{-30,0},{30,0}},
              color={158,66,200},
              thickness=0.5)},
          coordinateSystem(preserveAspectRatio=false)),
         Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>This connector can be used to connect a unidirectional outlet to a undirected fore port. </u>
<u>The state from the inlet is given to the backward direction of the rear port, the total potential as well as the massflow of inlet and port are set equal. </u>
<u>Note that when the flow is reversed, the resulting inertial potential can be different on both sides of this connector. </u>
</html>"));
    end ConnectInletRear;

    model ConnectRearOutlet
      "Directed/undirected connector with rear and outlet"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));
      parameter Boolean useDefaultStateAsRear = false "Use Default Medium states as state_rearwards";

      Chemical.Undirected.Interfaces.Rear rear(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
      Chemical.Interfaces.Outlet outlet(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={30,0}),
            iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));
      Chemical.Interfaces.StateInput state_rear if not useDefaultStateAsRear
        annotation (Placement(
            transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,40})));


    protected
      outer Chemical.DropOfCommons dropOfCommons;
       Modelica.Units.SI.ChemicalPotential r_rear, r_outlet "Inertial potential";

    equation

      connect(rear.state_rearwards, state_rear);

      outlet.state = rear.state_forwards;

      rear.n_flow + outlet.n_flow = 0;
      r_outlet + outlet.state.u = r_rear + Chemical.Undirected.Internal.regStep(
            rear.n_flow,
            rear.state_forwards.u,
            rear.state_rearwards.u,
            dropOfCommons.n_flow_reg);

      L/2*der(outlet.n_flow) = outlet.r - r_outlet;
      L/2*der(rear.n_flow) = rear.r - r_rear;

      if useDefaultStateAsRear then
        rear.state_rearwards = Chemical.Interfaces.SubstanceState(0,0);
      end if;


      outlet.solution = rear.solution;
      outlet.definition = rear.definition;
      annotation (Icon(
          graphics={
            Line(
              points={{-30,0},{30,0}},
              color={158,66,200},
              thickness=0.5),
            Line(points={{2,58},{0,58}}, color={158,66,200}),
            Line(
              points={{0,0},{0,60}},
              color={162,29,33},
              arrow={Arrow.Filled,Arrow.None},
              arrowSize = 20)},
          coordinateSystem(preserveAspectRatio=false)),
         Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>This connector can be used to connect a unidirectional inlet to a undirected fore port. </u>
<u>The state from the forward direction of the rear port is handed to the outlet, the total potential as well as the massflow of outlet and port are set equal. </u>
<u>Note that when the flow is reversed, the resulting inertial potential can be different on both sides of this connector. </u>
<u>Since the unidirectional side of this connector is an oulet, there is no information about the rear state from this side. This information must be given by the StateInlet.</u>
</html>"));
    end ConnectRearOutlet;

    model ConnectForeOutlet
      "Directed/undirected connector with rear and outlet"
      extends
        Chemical.Undirected.Topology.Internal.PartialSubstanceAndSolutionSource;

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));
      parameter Boolean useDefaultStateAsRear = false "Use Default Medium states as state_rearwards";

      Chemical.Undirected.Interfaces.Fore fore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-40,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
      Chemical.Interfaces.Outlet outlet(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={40,0}),
            iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));
      ConnectRearOutlet connectRearOutlet(redeclare package stateOfMatter =
            stateOfMatter, final L=L/2, final useDefaultStateAsRear = useDefaultStateAsRear)
        annotation (Placement(transformation(extent={{0,-10},{20,10}})));
      ConnectForeFore connectForeFore(redeclare package stateOfMatter =
            stateOfMatter,
        substanceData=substanceData, solutionParam=solutionParam, useSolution=useSolution,
                           final L=L/2)
        annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
      Chemical.Interfaces.StateInput state_rear if not useDefaultStateAsRear
        annotation (Placement(
            transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,40})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      if useSolution then
        connect(solution,connectForeFore.solution);
      end if;
      connect(connectRearOutlet.rear, connectForeFore.fore_b) annotation (Line(
          points={{7,0},{-7,0}},
          color={158,66,200},
          thickness=0.5));
      connect(connectRearOutlet.outlet, outlet) annotation (Line(
          points={{13,0},{40,0}},
          color={158,66,200},
          thickness=0.5));
      connect(connectForeFore.fore_a, fore) annotation (Line(
          points={{-13,0},{-40,0}},
          color={158,66,200},
          thickness=0.5));
      connect(connectRearOutlet.state_rear, state_rear) annotation (Line(points={{10,4},{
              10,12},{0,12},{0,40}}, color={162,29,33}));
      annotation (Icon(
          graphics={
            Line(
              points={{-30,0},{30,0}},
              color={158,66,200},
              thickness=0.5), Line(
              points={{0,0},{0,60}},
              color={162,29,33},
              arrow={Arrow.Filled,Arrow.None},
              arrowSize = 20)},
          coordinateSystem(preserveAspectRatio=false)),
         Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>This connector can be used to connect a unidirectional inlet to a undirected rear port. </u>
<u>The state from the rearward direction of the fore port is handed to the outlet, the total potential as well as the massflow of outlet and port are set equal. </u>
<u>Note that when the flow is reversed, the resulting inertial potential can be different on both sides of this connector. </u>
<u>Since the unidirectional side of this connector is an oulet, there is no information about the rear state from this side. This information must be given by the StateInlet.</u>
</html>"));
    end ConnectForeOutlet;

    model ConnectForeFore "Undirected connector with fore and fore"
      extends Chemical.Undirected.Topology.Internal.PartialSubstanceAndSolutionSource;
     /* replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby 
    Chemical.Interfaces.StateOfMatter
  "Substance model to translate data into substance properties"
    annotation (choices(
      choice(redeclare package stateOfMatter =
        Chemical.Interfaces.Incompressible  "Incompressible"),
      choice(redeclare package stateOfMatter =
        Chemical.Interfaces.IdealGas        "Ideal Gas"),
      choice(redeclare package stateOfMatter =
        Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
      choice(redeclare package stateOfMatter =
        Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
*/
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Fore fore_a(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
      Chemical.Undirected.Interfaces.Fore fore_b(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation
      fore_b.state_forwards = fore_a.state_rearwards;
      fore_a.state_forwards = fore_b.state_rearwards;
      fore_a.n_flow + fore_b.n_flow = 0;
      L*der(fore_a.n_flow) = fore_a.r - fore_b.r;

      fore_b.definition=substanceData;
      fore_a.definition=substanceData;
      fore_b.solution=solutionState;
      fore_a.solution=solutionState;

      annotation (Icon(
          graphics={
            Line(
              points={{-20,0},{20,0}},
              color={158,66,200},
              thickness=0.5)},
          coordinateSystem(preserveAspectRatio=false)),
         Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>This connector can be used to connect two undirected rear ports. </u>
<u>Basically the connector switches the names of output and input of the two ports.</u>
</html>"));
    end ConnectForeFore;

    model ConnectRearRear "Undirected connector with rear and rear"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Rear rear_a(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
      Chemical.Undirected.Interfaces.Rear rear_b(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation
      rear_b.state_rearwards = rear_a.state_forwards;
      rear_a.state_rearwards = rear_b.state_forwards;
      rear_a.n_flow + rear_b.n_flow = 0;
      L*der(rear_a.n_flow) = rear_a.r - rear_b.r;

      annotation (Icon(
          graphics={
            Line(
              points={{-20,0},{20,0}},
              color={158,66,200},
              thickness=0.5)},
          coordinateSystem(preserveAspectRatio=false)),
         Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>This connector can be used to connect two undirected fore ports. </u>
<u>Basically the connector switches the names of output and input of the two ports.</u>
</html>"));
    end ConnectRearRear;

    model JunctionRFF "Junction with rear and two fores"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Rear rear(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
      Chemical.Undirected.Interfaces.Fore foreA(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
      Chemical.Undirected.Interfaces.Fore foreB(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
      JunctionMN junctionMN(
        final M=2,
        final N=1,
        final n_flow_reg=n_flow_reg,
        final L=L,
        redeclare package stateOfMatter=stateOfMatter)
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(junctionMN.rears[1], rear) annotation (Line(
          points={{-50,0},{-76,0},{-76,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.fores[1], foreB) annotation (Line(
          points={{-30,-0.5},{0,-0.5},{0,-100}},
          color={158,66,200},
          thickness=0.5));
      connect(foreA, junctionMN.fores[2]) annotation (Line(
          points={{0,100},{0,0.5},{-30,0.5}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,-100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{20,100},{60,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{20,-60},{60,-100}},
              textColor={175,175,175},
              textString="B")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>
Junction with a rear and two fores in a lying T shape.
</u>
</html>"));
    end JunctionRFF;

    model JunctionRFF2 "Junction with rear and two fores"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Rear rear(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
      Chemical.Undirected.Interfaces.Fore foreA(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
      Chemical.Undirected.Interfaces.Fore foreB(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
      JunctionMN junctionMN(
        final M=2,
        final N=1,
        final n_flow_reg=n_flow_reg,
        final L=L,
        redeclare package stateOfMatter=stateOfMatter)
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(junctionMN.rears[1], rear) annotation (Line(
          points={{-50,0},{-76,0},{-76,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.fores[2], foreA) annotation (Line(
          points={{-30,0.5},{-30,1},{0,1},{0,100}},
          color={158,66,200},
          thickness=0.5));
      connect(foreB, junctionMN.fores[1]) annotation (Line(
          points={{100,0},{20,0},{20,-0.5},{-30,-0.5}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{20,100},{60,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{60,-20},{100,-60}},
              textColor={175,175,175},
              textString="B")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Junction with a rear and two fores in a standing T shape.</u>
</html>"));
    end JunctionRFF2;

    model JunctionRRF "Junction with two rears and a fore"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
       parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Fore fore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
      Chemical.Undirected.Interfaces.Rear rearA(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
      Chemical.Undirected.Interfaces.Rear rearB(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
      JunctionMN junctionMN(
        final M=1,
        final N=2,
        final n_flow_reg=n_flow_reg,
        final L=L,
        redeclare package stateOfMatter=stateOfMatter)
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(rearB, junctionMN.rears[1]) annotation (Line(
          points={{0,-100},{0,-2},{30,-2},{30,-0.5}},
          color={158,66,200},
          thickness=0.5));
      connect(rearA, junctionMN.rears[2]) annotation (Line(
          points={{0,100},{0,0.5},{30,0.5}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.fores[1], fore) annotation (Line(
          points={{50,0},{76,0},{76,0},{100,0}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,-100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{20,100},{60,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{20,-60},{60,-100}},
              textColor={175,175,175},
              textString="B")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Junction with two rears and a fore in a lying T shape.</u>
</html>"));
    end JunctionRRF;

    model JunctionRRF2 "Junction with two rears and a fore"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
       parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

      JunctionMN junctionMN(
        final M=1,
        final N=2,
        final n_flow_reg=n_flow_reg,
        final L=L,
        redeclare package stateOfMatter=stateOfMatter)
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));
      Chemical.Undirected.Interfaces.Rear rearA(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
      Chemical.Undirected.Interfaces.Rear rearB(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
      Chemical.Undirected.Interfaces.Fore fore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(junctionMN.rears[1], rearA) annotation (Line(
          points={{30,-0.5},{0,-0.5},{0,-100}},
          color={158,66,200},
          thickness=0.5));
      connect(fore, junctionMN.fores[1]) annotation (Line(
          points={{100,0},{76,0},{76,0},{50,0}},
          color={158,66,200},
          thickness=0.5));
      connect(rearB, junctionMN.rears[2]) annotation (Line(
          points={{-100,0},{-36,0},{-36,0.5},{30,0.5}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,-100}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{20,-100},{60,-60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{60,20},{100,60}},
              textColor={175,175,175},
              textString="B")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Junction with two rears and a fore in a standing T shape.</u>
</html>"));
    end JunctionRRF2;

    model JunctionRFFF "Junction with a rear and three fores"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Rear rear(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
      Chemical.Undirected.Interfaces.Fore foreA(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
      Chemical.Undirected.Interfaces.Fore foreB(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
      Chemical.Undirected.Interfaces.Fore foreC(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
      JunctionMN junctionMN(
        final M=3,
        final N=1,
        final n_flow_reg=n_flow_reg,
        final L=L,
        redeclare package stateOfMatter=stateOfMatter)
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation
      connect(junctionMN.rears[1], rear) annotation (Line(
          points={{-50,0},{-76,0},{-76,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.fores[3], foreA) annotation (Line(
          points={{-30,0.666667},{0,0.666667},{0,100}},
          color={158,66,200},
          thickness=0.5));
      connect(foreB, junctionMN.fores[2]) annotation (Line(
          points={{100,0},{-30,0}},
          color={158,66,200},
          thickness=0.5));
      connect(foreC, junctionMN.fores[1]) annotation (Line(
          points={{0,-100},{0,-0.666667},{-30,-0.666667}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,-100},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{-20,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{50,20},{90,60}},
              textColor={175,175,175},
              textString="B"),
            Text(
              extent={{60,-100},{20,-60}},
              textColor={175,175,175},
              textString="C")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Junction with a rear and three fores in a x shape.</u>
</html>"));
    end JunctionRFFF;

    model JunctionRRFF "Junction with two rears and two fores"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Rear reara(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
      Chemical.Undirected.Interfaces.Rear rearb(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
      Chemical.Undirected.Interfaces.Fore foreB(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
      Chemical.Undirected.Interfaces.Fore foreA(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
      JunctionMN junctionMN(
        final M=2,
        final N=2,
        final n_flow_reg=n_flow_reg,
        final L=L,
        redeclare package stateOfMatter=stateOfMatter)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation
      connect(foreB, junctionMN.fores[2]) annotation (Line(
          points={{100,0},{40,0},{40,0.5},{10,0.5}},
          color={158,66,200},
          thickness=0.5));
      connect(foreA, junctionMN.fores[1]) annotation (Line(
          points={{0,-100},{0,-20},{10,-20},{10,-0.5}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.rears[1], reara) annotation (Line(
          points={{-10,-0.5},{-40,-0.5},{-40,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(rearb, junctionMN.rears[2]) annotation (Line(
          points={{0,100},{0,20},{-10,20},{-10,0.5}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,-100},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{-20,60}},
              textColor={175,175,175},
              textString="b"),
            Text(
              extent={{50,20},{90,60}},
              textColor={175,175,175},
              textString="B"),
            Text(
              extent={{60,-100},{20,-60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{-50,-20},{-90,-60}},
              textColor={175,175,175},
              textString="a")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Junction with two rears and two fores in a x shape.</u>
</html>"));
    end JunctionRRFF;

    model JunctionRRFF2 "Junction with two rears and two fores"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Rear reara(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
      Chemical.Undirected.Interfaces.Rear rearb(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
      Chemical.Undirected.Interfaces.Fore foreA(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
      Chemical.Undirected.Interfaces.Fore foreB(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
      JunctionMN junctionMN(
        final M=2,
        final N=2,
        final n_flow_reg=n_flow_reg,
        final L=L,
        redeclare package stateOfMatter=stateOfMatter)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}},
            rotation=90,
            origin={0,10})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation
      connect(rearb, junctionMN.rears[1]) annotation (Line(
          points={{0,-100},{0,-20},{20,-20},{20,20},{0.5,20}},
          color={158,66,200},
          thickness=0.5));
      connect(reara, junctionMN.rears[2]) annotation (Line(
          points={{0,100},{0,20},{-0.5,20}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.fores[2], foreA) annotation (Line(
          points={{-0.5,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(foreB, junctionMN.fores[1]) annotation (Line(
          points={{100,0},{0.5,0}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,-100},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{-20,60}},
              textColor={175,175,175},
              textString="a"),
            Text(
              extent={{50,20},{90,60}},
              textColor={175,175,175},
              textString="B"),
            Text(
              extent={{60,-100},{20,-60}},
              textColor={175,175,175},
              textString="b"),
            Text(
              extent={{-50,-20},{-90,-60}},
              textColor={175,175,175},
              textString="A")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Junction with two rears and two fores in a x shape.</u>
</html>"));
    end JunctionRRFF2;

    model JunctionRRRF "Junction with tree rears and a fore"

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Rear rearA(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
      Chemical.Undirected.Interfaces.Fore fore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
      JunctionMN junctionMN(
        final M=1,
        final N=3,
        final n_flow_reg=n_flow_reg,
        final L=L,
        redeclare package stateOfMatter=stateOfMatter)
        annotation (Placement(transformation(extent={{0,-10},{20,10}})));
      Chemical.Undirected.Interfaces.Rear rearB(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
      Chemical.Undirected.Interfaces.Rear rearC(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation
      connect(rearA, junctionMN.rears[2]) annotation (Line(
          points={{0,100},{0,0}},
          color={158,66,200},
          thickness=0.5));
      connect(rearB, junctionMN.rears[1]) annotation (Line(
          points={{-100,0},{0,0},{0,-0.666667}},
          color={158,66,200},
          thickness=0.5));
      connect(rearC, junctionMN.rears[3]) annotation (Line(
          points={{0,-100},{0,0.666667}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.fores[1], fore) annotation (Line(
          points={{20,0},{60,0},{60,0},{100,0}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,-100},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{-20,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{-50,-20},{-90,-60}},
              textColor={175,175,175},
              textString="B"),
            Text(
              extent={{60,-100},{20,-60}},
              textColor={175,175,175},
              textString="C")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Junction with three rears and a fore in a x shape.</u>
</html>"));
    end JunctionRRRF;

    model ConnectorInletOutletFore

     replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_ref=dropOfCommons.n_flow_reg "Reference mass flow" annotation (Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.ChemicalPotential u_ref=1e5 "Reference potential" annotation (Dialog(tab="Advanced"));
      parameter Boolean assumeConstantDensity = true "If true only mass-flow rate will determine the mixing"
        annotation (Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));

      Chemical.Undirected.Interfaces.Fore fore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-120},{20,-80}}), iconTransformation(extent={{-20,-120},{20,-80}})));
      Chemical.Interfaces.Inlet inlet(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{
                -80,20}})));
      Chemical.Interfaces.Outlet outlet(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,
                20}})));
      Chemical.FlowControl.CheckValve checkValve(
        redeclare package stateOfMatter = stateOfMatter,
        final L=L,
        final n_flow_ref=n_flow_ref,
        final u_ref=u_ref)
          annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
      Chemical.FlowControl.CheckValve checkValve1(
        redeclare package stateOfMatter = stateOfMatter,
        final L=L,
        final n_flow_ref=n_flow_ref,
        final u_ref=u_ref)
          annotation (Placement(transformation(extent={{50,-10},{70,10}})));
      ConnectRearOutlet connectRearOutlet(
        redeclare package stateOfMatter = stateOfMatter,
        final L=L,
        final useDefaultStateAsRear=false)
          annotation (Placement(transformation(extent={{20,-10},{40,10}})));
      ConnectInletFore connectInletFore(
        redeclare package stateOfMatter = stateOfMatter,
        final L=L)
          annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
      JunctionRFF2 junctionRFF2_1(
        redeclare package stateOfMatter = stateOfMatter,
        final assumeConstantDensity=assumeConstantDensity,
        final n_flow_reg=n_flow_reg,
        final L=L)
          annotation (Placement(transformation(extent={{-10,10},{10,-10}})));
      Chemical.Sensors.SensorState sensorState(redeclare package stateOfMatter =
            stateOfMatter)
          annotation (Placement(transformation(extent={{-10,8},{10,28}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation
      connect(junctionRFF2_1.foreA, fore) annotation (Line(
          points={{0,-10},{0,-52},{0,-100},{0,-100}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionRFF2_1.foreB, connectRearOutlet.rear)
        annotation (Line(
          points={{10,0},{27,0}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionRFF2_1.rear, connectInletFore.fore) annotation (Line(
          points={{-10,0},{-27,0}},
          color={158,66,200},
          thickness=0.5));
      connect(connectInletFore.inlet, checkValve.outlet) annotation (Line(
          points={{-33,0},{-50,0}},
          color={158,66,200},
          thickness=0.5));
      connect(checkValve.inlet, inlet) annotation (Line(
          points={{-70,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(connectRearOutlet.outlet, checkValve1.inlet) annotation (Line(
          points={{33,0},{50,0}},
          color={158,66,200},
          thickness=0.5));
      connect(checkValve1.outlet, outlet) annotation (Line(
          points={{70,0},{100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(sensorState.inlet, checkValve.outlet)
        annotation (Line(
          points={{-10,18},{-40,18},{-40,0},{-50,0}},
          color={158,66,200},
          thickness=0.5));
      connect(sensorState.state_out, connectRearOutlet.state_rear) annotation (Line(points={{10,18},{30,18},{30,4}}, color={162,29,33}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,-100},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5)}), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end ConnectorInletOutletFore;

    package Tests "Test models for the undirected topology package"
    extends Modelica.Icons.ExamplesPackage;

      model TestJunction "Test for the undirected junction"
        extends Modelica.Icons.Example;

        replaceable package Medium = Chemical.Media.myMedia.Air.SimpleAir constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medum model for the Test" annotation (Documentation(info="<html>
<u>This is the replaceable package that determines the medium of the Test. </u>
</html>"));

        replaceable function pLoss =
            Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
            (
          material=Chemical.Processes.Internal.Material.wood)
          constrainedby Chemical.Processes.Internal.FlowResistance.partialChemicalPotentialLoss
          "ChemicalPotential loss function for all Flow resistances"
          annotation(choicesAllMatching = true, Documentation(info="<html>
<u>This is the potential loss function used for all resistances except the two on the outlets of the right two cases.</u>
</html>"));

        Chemical.Boundaries.Source source(redeclare package Medium=Medium,
            u0_par=200000)
          annotation (Placement(transformation(extent={{-160,30},{-140,50}})));
        Chemical.Boundaries.Sink sink(redeclare package Medium=Medium,
            u0_par=150000)
          annotation (Placement(transformation(extent={{-20,50},{0,70}})));
        Chemical.Boundaries.Sink sink1(redeclare package Medium=Medium,
            u0_par=100000)
          annotation (Placement(transformation(extent={{-20,10},{0,30}})));
        Chemical.Topology.SplitterT1 splitterT1_1(redeclare package Medium=Medium, L=
              2*dropOfCommons.L)
          annotation (Placement(transformation(extent={{-94,30},{-74,50}})));
        inner Chemical.DropOfCommons dropOfCommons(L=1e2) annotation (Placement(transformation(extent={{-134,-8},{-114,12}})));
        Chemical.Processes.FlowResistance flowResistance(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{-46,50},{-26,70}})));
        Chemical.Processes.FlowResistance flowResistance1(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{-46,10},{-26,30}})));
        Chemical.Boundaries.Source source2(
          redeclare package stateOfMatter = stateOfMatter, u0_par=2000000)
          annotation (Placement(transformation(extent={{0,50},{20,70}})));
        Chemical.Boundaries.Source source3(
          redeclare package stateOfMatter = stateOfMatter, u0_par=3500000)
          annotation (Placement(transformation(extent={{0,10},{20,30}})));
        Chemical.Processes.FlowResistance flowResistance4(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{26,10},{46,30}})));
        Chemical.Processes.FlowResistance flowResistance5(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{26,50},{46,70}})));
        Chemical.Topology.JunctionT1 junctionT1_1(redeclare package Medium=Medium,
          assumeConstantDensity=false,
          L=2*dropOfCommons.L)
          annotation (Placement(transformation(extent={{94,30},{74,50}})));
        Chemical.Boundaries.Sink sink4(
          redeclare package stateOfMatter = stateOfMatter, u0_par=100000)
          annotation (Placement(transformation(extent={{140,30},{160,50}})));
        Chemical.Boundaries.Source source6(redeclare package stateOfMatter =
              stateOfMatter,
            u0_par=200000)
          annotation (Placement(transformation(extent={{-160,-50},{-140,-30}})));
        Chemical.Boundaries.Sink sink6(redeclare package stateOfMatter =
              stateOfMatter,
            u0_par=150000)
          annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
        Chemical.Boundaries.Sink sink7(redeclare package stateOfMatter =
              stateOfMatter,
            u0_par=100000)
          annotation (Placement(transformation(extent={{-20,-70},{0,-50}})));
        ConnectInletFore connectInletFore1(redeclare package stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-118,-50},{-98,-30}})));
        ConnectRearOutlet connectRearOutlet2(redeclare package stateOfMatter =
              stateOfMatter,
            useDefaultStateAsRear=true)
          annotation (Placement(transformation(extent={{-66,-30},{-46,-10}})));
        ConnectRearOutlet connectRearOutlet3(redeclare package stateOfMatter =
              stateOfMatter,
            useDefaultStateAsRear=true)
          annotation (Placement(transformation(extent={{-66,-70},{-46,-50}})));
        Chemical.Processes.FlowResistance flowResistance8(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{-46,-30},{-26,-10}})));
        Chemical.Processes.FlowResistance flowResistance9(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{-46,-70},{-26,-50}})));
        Chemical.Boundaries.Source source7(redeclare package stateOfMatter =
              stateOfMatter,                                                                u0_par=2000000)
          annotation (Placement(transformation(extent={{0,-30},{20,-10}})));
        Chemical.Boundaries.Source source8(redeclare package stateOfMatter =
              stateOfMatter,                                                                u0_par=3500000)
          annotation (Placement(transformation(extent={{0,-70},{20,-50}})));
        Chemical.Processes.FlowResistance flowResistance10(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{26,-70},{46,-50}})));
        Chemical.Processes.FlowResistance flowResistance11(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{26,-30},{46,-10}})));
        Chemical.Boundaries.Sink sink8(redeclare package stateOfMatter =
              stateOfMatter,
            u0_par=100000)
          annotation (Placement(transformation(extent={{140,-50},{160,-30}})));
        ConnectInletRear connectInletRear2(redeclare package stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{46,-30},{66,-10}})));
        ConnectInletRear connectInletRear3(redeclare package stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{46,-70},{66,-50}})));
        ConnectForeOutlet connectForeOutlet1(redeclare package stateOfMatter =
              stateOfMatter,
            useDefaultStateAsRear=true)
          annotation (Placement(transformation(extent={{98,-50},{118,-30}})));
        JunctionMN junctionMN(redeclare package stateOfMatter = stateOfMatter,
          M=2,
          N=1)
          annotation (Placement(transformation(extent={{-94,-50},{-74,-30}})));
        JunctionMN junctionMN1(redeclare package stateOfMatter = stateOfMatter,
          M=2,
          N=1,
          assumeConstantDensity=false)
          annotation (Placement(transformation(extent={{94,-50},{74,-30}})));
        Chemical.Processes.FlowResistance flowResistance2(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{-136,30},{-116,50}})));
        Chemical.Processes.FlowResistance flowResistance3(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss = pLoss)
          annotation (Placement(transformation(extent={{-136,-50},{-116,-30}})));
        Chemical.Processes.FlowResistance flowResistance6(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (k=1000))
          annotation (Placement(transformation(extent={{116,-50},{136,-30}})));
        Chemical.Processes.FlowResistance flowResistance7(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.1,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (k=1000))
          annotation (Placement(transformation(extent={{116,30},{136,50}})));
      equation
        connect(sink.inlet, flowResistance.outlet) annotation (Line(
            points={{-20,60},{-26,60}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance.inlet, splitterT1_1.outletA) annotation (Line(
            points={{-46,60},{-84,60},{-84,50}},
            color={158,66,200},
            thickness=0.5));
        connect(sink1.inlet, flowResistance1.outlet) annotation (Line(
            points={{-20,20},{-26,20}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance1.inlet, splitterT1_1.outletB) annotation (Line(
            points={{-46,20},{-84,20},{-84,30}},
            color={158,66,200},
            thickness=0.5));
        connect(source2.outlet, flowResistance5.inlet) annotation (Line(
            points={{20,60},{26,60}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance4.inlet, source3.outlet) annotation (Line(
            points={{26,20},{20,20}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance4.outlet, junctionT1_1.inletB) annotation (Line(
            points={{46,20},{84,20},{84,30}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionT1_1.inletA, flowResistance5.outlet) annotation (Line(
            points={{84,50},{84,60},{46,60}},
            color={158,66,200},
            thickness=0.5));
        connect(sink6.inlet,flowResistance8. outlet) annotation (Line(
            points={{-20,-20},{-26,-20}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance8.inlet, connectRearOutlet2.outlet) annotation (Line(
            points={{-46,-20},{-53,-20}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance9.inlet,connectRearOutlet3. outlet) annotation (Line(
            points={{-46,-60},{-53,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(sink7.inlet,flowResistance9. outlet) annotation (Line(
            points={{-20,-60},{-26,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance10.inlet, source8.outlet) annotation (Line(
            points={{26,-60},{20,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(source7.outlet, flowResistance11.inlet) annotation (Line(
            points={{20,-20},{26,-20}},
            color={158,66,200},
            thickness=0.5));
        connect(connectInletRear2.inlet, flowResistance11.outlet) annotation (Line(
            points={{53,-20},{46,-20}},
            color={158,66,200},
            thickness=0.5));
        connect(connectInletRear3.inlet, flowResistance10.outlet) annotation (Line(
            points={{53,-60},{46,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(connectRearOutlet2.rear, junctionMN.fores[2]) annotation (Line(
            points={{-59,-20},{-66,-20},{-66,-39.5},{-74,-39.5}},
            color={158,66,200},
            thickness=0.5));
        connect(connectRearOutlet3.rear, junctionMN.fores[1]) annotation (Line(
            points={{-59,-60},{-66,-60},{-66,-40.5},{-74,-40.5}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionMN.rears[1], connectInletFore1.fore) annotation (Line(
            points={{-94,-40},{-100,-40},{-100,-40},{-105,-40}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionMN1.rears[1], connectForeOutlet1.fore) annotation (Line(
            points={{94,-40},{100,-40},{100,-40},{105,-40}},
            color={158,66,200},
            thickness=0.5));
        connect(connectInletRear2.rear, junctionMN1.fores[2]) annotation (Line(
            points={{59,-20},{66,-20},{66,-39.5},{74,-39.5}},
            color={158,66,200},
            thickness=0.5));
        connect(connectInletRear3.rear, junctionMN1.fores[1]) annotation (Line(
            points={{59,-60},{66,-60},{66,-40.5},{74,-40.5}},
            color={158,66,200},
            thickness=0.5));
        connect(source.outlet, flowResistance2.inlet) annotation (Line(
            points={{-140,40},{-136,40}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT1_1.inlet, flowResistance2.outlet) annotation (Line(
            points={{-94,40},{-116,40}},
            color={158,66,200},
            thickness=0.5));
        connect(connectInletFore1.inlet, flowResistance3.outlet) annotation (Line(
            points={{-111,-40},{-116,-40}},
            color={158,66,200},
            thickness=0.5));
        connect(source6.outlet, flowResistance3.inlet) annotation (Line(
            points={{-140,-40},{-136,-40}},
            color={158,66,200},
            thickness=0.5));
        connect(connectForeOutlet1.outlet, flowResistance6.inlet) annotation (Line(
            points={{111,-40},{116,-40}},
            color={158,66,200},
            thickness=0.5));
        connect(sink8.inlet, flowResistance6.outlet) annotation (Line(
            points={{140,-40},{136,-40}},
            color={158,66,200},
            thickness=0.5));
        connect(sink4.inlet, flowResistance7.outlet) annotation (Line(
            points={{140,40},{136,40}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionT1_1.outlet, flowResistance7.inlet) annotation (Line(
            points={{94,40},{116,40}},
            color={158,66,200},
            thickness=0.5));
        annotation (
        experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
        Documentation(info="<html>
<u>This model tests the undirected junction against the unidirectional junction. The states of the left two and the right two systems are expected to be the same, when the junctions have the same settings. </u>
<u>Note that the unidirectional junctions have two times the inertance, since the undirected junction comes with a additional connector, which in turn adds inertance to each outlet.</u>
<u><br>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"),Diagram(coordinateSystem(extent={{-160,-80},{160,80}})),
          Icon(coordinateSystem(extent={{-100,-100},{100,100}})));
      end TestJunction;

      model TestConnectors "Test for the connectors"
        extends Modelica.Icons.Example;

        replaceable package Medium = Chemical.Media.myMedia.Air.SimpleAir constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
                                                               "Medum model for the Test"
          annotation (Documentation(info="<html>
<u>This is the replaceable package that determines the medium of the Test. </u>
</html>"));

        ConnectRearRear connectRearRear(redeclare package Medium=Medium)
          annotation (Placement(transformation(extent={{-10,60},{10,80}})));
        ConnectForeFore connectForeFore(redeclare package Medium=Medium)
          annotation (Placement(transformation(extent={{-10,40},{10,60}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package stateOfMatter = stateOfMatter,
          u0_par=100000,
          fore(n_flow(start=0, fixed=true))) annotation (Placement(transformation(extent={{-40,60},{-20,80}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=100000)
          annotation (Placement(transformation(extent={{-20,40},{-40,60}})));
        Chemical.Boundaries.Source source(redeclare package Medium=Medium,
            u0_par=100000)
          annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
        ConnectInletFore connectInletFore(redeclare package stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-10,20},{10,40}})));
        ConnectRearOutlet connectRearOutlet(redeclare package Medium=Medium, rear(n_flow(start=0, fixed=true)))
          annotation (Placement(transformation(extent={{-10,0},{10,20}})));
        Chemical.Boundaries.Sink sink(redeclare package Medium=Medium,
            potentialFromInput=true)
          annotation (Placement(transformation(extent={{20,0},{40,20}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear1(redeclare
            package                                                                  stateOfMatter =
              stateOfMatter,                                                                                        potentialFromInput=true)
          annotation (Placement(transformation(extent={{40,60},{20,80}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore1(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          rear(n_flow(start=0, fixed=true))) annotation (Placement(transformation(extent={{20,40},{40,60}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore2(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          rear(n_flow(start=0, fixed=true))) annotation (Placement(transformation(extent={{20,20},{40,40}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear2(redeclare
            package                                                                  stateOfMatter =
              stateOfMatter,                                                                                        u0_par=100000)
          annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
        Chemical.Boundaries.Source source1(redeclare package stateOfMatter =
              stateOfMatter,
          u0_par=100000,
          L=1.5*dropOfCommons.L,
          outlet(n_flow(start=0, fixed=true)))
          annotation (Placement(transformation(extent={{-40,100},{-20,120}})));
        Chemical.Boundaries.Sink sink1(redeclare package stateOfMatter =
              stateOfMatter,
          potentialFromInput=true,
          L=1.5*dropOfCommons.L)
          annotation (Placement(transformation(extent={{20,100},{40,120}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear3(
          redeclare package stateOfMatter = stateOfMatter,
          u0_par=100000,
          L=1.5*dropOfCommons.L) annotation (Placement(transformation(extent={{-40,80},{-20,100}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore3(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          L=1.5*dropOfCommons.L,
          rear(n_flow(start=0, fixed=true))) annotation (Placement(transformation(extent={{20,80},{40,100}})));
        inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{-88,80},{-68,100}})));
        Modelica.Blocks.Sources.Pulse pulse(
          amplitude=1e5,
          period=5,
          offset=0.5e5,
          startTime=1)
          annotation (Placement(transformation(extent={{82,50},{62,70}})));
        Chemical.Boundaries.Source source2(redeclare package stateOfMatter =
              stateOfMatter,
            u0_par=100000)
          annotation (Placement(transformation(extent={{-50,-40},{-30,-20}})));
        ConnectInletFore connectInletFore1(redeclare package stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-30,-40},{-10,-20}})));
        ConnectRearOutlet connectRearOutlet1(redeclare package stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{10,-40},{30,-20}})));
        Chemical.Boundaries.Sink sink2(redeclare package stateOfMatter =
              stateOfMatter,
            potentialFromInput=true)
          annotation (Placement(transformation(extent={{30,-40},{50,-20}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          computeL=false,
          r=0.01,
          l=1,
          redeclare function pLoss =
              .Chemical.Processes.Internal.FlowResistance.laminarChemicalPotentialLoss)
          annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance2(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          computeL=false,
          r=0.01,
          l=1,
          redeclare function pLoss =
              .Chemical.Processes.Internal.FlowResistance.laminarChemicalPotentialLoss)
          annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear4(redeclare
            package                                                                  stateOfMatter =
              stateOfMatter,                                                                                        u0_par=100000)
          annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore4(redeclare
            package                                                                  stateOfMatter =
              stateOfMatter,                                                                                        potentialFromInput=true)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,-50})));
        Chemical.Boundaries.CreateState createState(
          redeclare package stateOfMatter = stateOfMatter,
          PFromInput=true)
          annotation (Placement(transformation(extent={{38,-12},{32,-4}})));
        ConnectorInletOutletFore connectorInletOutletFore(redeclare package stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={0,-68})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore5(redeclare
            package                                                                  stateOfMatter =
              stateOfMatter,                                                                                        potentialFromInput=true)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={64,-90})));
        Chemical.Boundaries.Source source3(redeclare package stateOfMatter =
              stateOfMatter,                                                                u0_par=200000)
          annotation (Placement(transformation(extent={{-80,-78},{-60,-58}})));
        Chemical.Boundaries.Sink sink3(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=false,
          u0_par=100000)
          annotation (Placement(transformation(extent={{54,-78},{74,-58}})));
        Chemical.Processes.FlowResistance flowResistance1(
          redeclare package stateOfMatter = stateOfMatter,
          r(displayUnit="mm") = 0.005,
          l=5,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss)
          annotation (Placement(transformation(extent={{-40,-78},{-20,-58}})));
        Chemical.Processes.FlowResistance flowResistance3(
          redeclare package stateOfMatter = stateOfMatter, initM_flow = Chemical.Utilities.Types.InitializationMethods.state,
          l=5,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss,
          r(displayUnit="mm") = 0.005)
          annotation (Placement(transformation(extent={{20,-78},{40,-58}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance4(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          l=5,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss,
          r(displayUnit="mm") = 0.005) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={30,-90})));
        Modelica.Blocks.Sources.Ramp ramp(
          height=1.4e5,
          duration=1,
          offset=0.8e5,
          startTime=4.5) annotation (Placement(transformation(extent={{100,-100},{80,-80}})));
      equation
        connect(boundary_rear.fore, connectRearRear.rear_a) annotation (Line(
            points={{-20,70},{-3,70}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore.rear, connectForeFore.fore_a) annotation (Line(
            points={{-20,50},{-3,50}},
            color={158,66,200},
            thickness=0.5));
        connect(connectInletFore.inlet, source.outlet) annotation (Line(
            points={{-3,30},{-20,30}},
            color={158,66,200},
            thickness=0.5));
        connect(sink.inlet, connectRearOutlet.outlet) annotation (Line(
            points={{20,10},{3,10}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear1.fore, connectRearRear.rear_b) annotation (Line(
            points={{20,70},{3,70}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore1.rear, connectForeFore.fore_b) annotation (Line(
            points={{20,50},{3,50}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore2.rear, connectInletFore.fore) annotation (Line(
            points={{20,30},{3,30}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear2.fore, connectRearOutlet.rear) annotation (Line(
            points={{-20,10},{-3,10}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear3.fore, boundary_fore3.rear) annotation (Line(
            points={{-20,90},{20,90}},
            color={158,66,200},
            thickness=0.5));
        connect(sink1.inlet, source1.outlet) annotation (Line(
            points={{20,110},{-20,110}},
            color={158,66,200},
            thickness=0.5));
        connect(pulse.y, sink1.u0_var) annotation (Line(points={{61,60},{46,60},{46,110},{32,110}},
                                                                                                color={0,0,127}));
        connect(pulse.y, boundary_fore3.u0_var) annotation (Line(points={{61,60},{46,60},{46,96},{32,96}},
                                                                                                         color={0,0,127}));
        connect(pulse.y, boundary_rear1.u0_var) annotation (Line(points={{61,60},{46,60},{46,76},{32,76}},
                                                                                                       color={0,0,127}));
        connect(pulse.y, boundary_fore1.u0_var) annotation (Line(points={{61,60},{46,60},{46,56},{32,56}}, color={0,0,127}));
        connect(pulse.y, boundary_fore2.u0_var) annotation (Line(points={{61,60},{46,60},{46,36},{32,36}}, color={0,0,127}));
        connect(pulse.y, sink.u0_var) annotation (Line(points={{61,60},{46,60},{46,10},{32,10}}, color={0,0,127}));
        connect(connectInletFore1.inlet, source2.outlet) annotation (Line(
            points={{-23,-30},{-30,-30}},
            color={158,66,200},
            thickness=0.5));
        connect(sink2.inlet, connectRearOutlet1.outlet) annotation (Line(
            points={{30,-30},{23,-30}},
            color={158,66,200},
            thickness=0.5));
        connect(connectInletFore1.fore, flowResistance.rear) annotation (Line(
            points={{-17,-30},{-10,-30}},
            color={158,66,200},
            thickness=0.5));
        connect(connectRearOutlet1.rear, flowResistance.fore) annotation (Line(
            points={{17,-30},{10,-30}},
            color={158,66,200},
            thickness=0.5));
        connect(sink2.u0_var, sink.u0_var) annotation (Line(points={{42,-30},{46,-30},{46,10},{32,10}},
                                    color={0,0,127}));
        connect(flowResistance2.rear, boundary_rear4.fore) annotation (Line(
            points={{-10,-50},{-30,-50}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance2.fore, boundary_fore4.rear) annotation (Line(
            points={{10,-50},{30,-50}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore4.u0_var, sink.u0_var) annotation (Line(points={{42,-44},{46,-44},{46,10},{32,10}},
                                              color={0,0,127}));
        connect(createState.y, connectRearOutlet.state_rear) annotation (Line(points={{32,-8},{12,-8},{12,20},{0,20},{0,14}},
                                                                              color={
                162,29,33}));
        connect(createState.u_inp, sink.u0_var) annotation (Line(points={{38,-4},{46,-4},{46,10},{32,10}},
                                         color={0,0,127}));
        connect(connectRearOutlet1.state_rear, connectRearOutlet.state_rear)
          annotation (Line(points={{20,-26},{20,-8},{12,-8},{12,28},{0,28},{0,14}},
              color={162,29,33}));
        connect(source3.outlet, flowResistance1.inlet) annotation (Line(
            points={{-60,-68},{-40,-68}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance1.outlet, connectorInletOutletFore.inlet)
          annotation (Line(
            points={{-20,-68},{-10,-68}},
            color={158,66,200},
            thickness=0.5));
        connect(sink3.inlet, flowResistance3.outlet) annotation (Line(
            points={{54,-68},{40,-68}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance3.inlet, connectorInletOutletFore.outlet)
          annotation (Line(
            points={{20,-68},{10,-68}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore5.rear, flowResistance4.fore) annotation (Line(
            points={{54,-90},{40,-90}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance4.rear, connectorInletOutletFore.fore)
          annotation (Line(
            points={{20,-90},{0,-90},{0,-78}},
            color={158,66,200},
            thickness=0.5));
        connect(ramp.y, boundary_fore5.u0_var) annotation (Line(points={{79,-90},{72,-90},{72,-84},{66,-84}},
                                                                                            color={0,0,127}));
        annotation (
        experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
        Documentation(info="<html>
<u>This model tests the connectors against the case when source is directly connected to sink. All massflows are expected to be the same, when the sources and sinks have the same settings. </u>
<u>Note that the sources and sinks without a connector have 1.5 times the inertance, since the additional connector adds inertance to each other path.</u>
<u><br>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"),Diagram(coordinateSystem(extent={{-100,-140},{100,120}})),
          Icon(coordinateSystem(extent={{-100,-100},{100,100}})));
      end TestConnectors;
      annotation (Documentation(info="<html>
<u>This package contains the test cases for the undirected topology package. </u>
</html>"));
    end Tests;

    package Internal
      partial model PartialSubstanceAndSolutionSource
        "Substance properties for components, where the substance is connected with the solution"

       replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
          Chemical.Interfaces.StateOfMatter
        "Substance model to translate data into substance properties"
          annotation (choices(
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.Incompressible  "Incompressible"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGas        "Ideal Gas"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

        parameter stateOfMatter.SubstanceDataParameters substanceData
       "Definition of the substance"
          annotation (choicesAllMatching = true, Dialog(enable=not useRear));

        parameter Boolean useSolution = false "Use solution connector?"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

        parameter Chemical.Interfaces.SolutionStateParameters solutionParam "Constant chemical solution state if not from rear or input"
          annotation (Dialog(enable=not useSolution and not useRear));

        Chemical.Interfaces.SolutionPort solution(
            T=solutionState.T,
            p=solutionState.p,
            v=solutionState.v,
            n=solutionState.n,
            m=solutionState.m,
            V=solutionState.V,
            G=solutionState.G,
            Q=solutionState.Q,
            I=solutionState.I,
            i=0,
            dH=0,
            dV=0,
            nj=0,
            mj=0,
            Vj=0,
            Gj=0,
            Qj=0,
            Ij=0)
              if useSolution "To connect substance with solution, where is pressented"
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

           Chemical.Interfaces.SolutionState solutionState;
      equation

        if not useSolution then
          solutionState.T=solutionParam.T "Temperature of the solution";
          solutionState.p=solutionParam.p "Pressure of the solution";
          solutionState.v=solutionParam.v "Electric potential in the solution";
          solutionState.n=solutionParam.n "Amount of the solution";
          solutionState.m=solutionParam.m "Mass of the solution";
          solutionState.V=solutionParam.V "Volume of the solution";
          solutionState.G=solutionParam.G "Free Gibbs energy of the solution";
          solutionState.Q=solutionParam.Q "Electric charge of the solution";
          solutionState.I=solutionParam.I "Mole fraction based ionic strength of the solution";
        end if;

      end PartialSubstanceAndSolutionSource;
    end Internal;
  annotation (Documentation(info="<html>
<u>This package contains the undirected junctions and necessary connectors between two undirected components, as well as directed-undirected connectors.</u>
<u>Note that in the undirected case it a distinction between junction and splitter is not possible.</u>
</html>",   revisions="<html>
<u><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</u>
</html>"),   Icon(graphics={
          Line(
            points={{-80,2},{12,2}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{12,2},{80,-78}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{12,2},{80,82}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{6,8},{18,-4}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)}));
  end Topology;

  package Processes "Undirected process package"
    model Reaction "Chemical Reaction"
      extends Chemical.Undirected.Processes.Internal.PartialReactionWithSubstanceData;
      extends Chemical.Interfaces.ConditionalKinetics
                                            (k_forward=1);
      Real rr_fore_exact,rr_rear_exact,  kb;
    equation

      rr = kf * Sx_fore * ( 1  -  exp(-duRT_fore));
      if nS>0 then
        rr = -kb * Px_rear * ( 1  -  exp(duRT_rear));
      end if;

      Kx = kb/kf;

      //the same as:
      rr_fore_exact = (kf*Sx_fore - kb*Px_fore);
      rr_rear_exact = (kf*Sx_rear - kb*Px_rear);



      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Rectangle(
              extent={{-100,-30},{100,30}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-100,-72},{100,-40}},
              lineColor={128,0,255},
            textString="%name"),
            Polygon(
              points={{-60,6},{-60,4},{54,4},{54,4},{18,14},{18,6},{-60,6}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{54,-8},{54,-6},{-60,-6},{-60,-6},{-24,-16},{-24,-8},{54,-8}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid)}),
        Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
<p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
<td><p>dissociation constant</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>G = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>G<sub>i</sub>) = &Delta;<sub>r</sub>H - T&middot;&Delta;<sub>r</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K) </p></td>
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
<h4><span style=\"color:#008000\">Notations</span></h4>
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
    extends Modelica.Icons.Package;

    model ForwardReaction "Chemical Reaction"
      extends Internal.PartialReactionWithSubstanceData;
      extends Chemical.Interfaces.ConditionalKinetics
                                            (k_forward=1);

    equation

      rr = kf * Sx_fore;

      if nP>0 then
         der(products[1].state_forwards.u).*TC = products[1].r;
      end if;



      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Rectangle(
              extent={{-100,-30},{100,30}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-100,-72},{100,-40}},
              lineColor={128,0,255},
            textString="%name"),
            Polygon(
              points={{-60,2},{-60,0},{54,0},{54,0},{18,10},{18,2},{-60,2}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-60,-2},{-60,0},{54,0},{54,0},{18,-10},{18,-2},{-60,-2}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid)}),
        Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
<p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
<td><p>dissociation constant</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>G = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>G<sub>i</sub>) = &Delta;<sub>r</sub>H - T&middot;&Delta;<sub>r</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K) </p></td>
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
<h4><span style=\"color:#008000\">Notations</span></h4>
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
    end ForwardReaction;

    model FlowResistance "Flow resistance model"
      extends Chemical.Undirected.Interfaces.SISOBiFlow(
                                    final L=if computeL then l/(r^2*pi) else L_value, final cliu_u_out=true);

      import Modelica.Constants.pi "Constant Pi";

      parameter Modelica.Units.SI.Radius r(min=0) "Radius of pipe";
      parameter Modelica.Units.SI.Length l(min=0) "Length of pipe";
      parameter Chemical.Utilities.Units.Inertance L_value=dropOfCommons.L "Inertance of pipe" annotation (Dialog(enable=not computeL, tab="Advanced"));
      parameter Boolean computeL = true "Compute L from r and l"
        annotation(Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimal input density" annotation (Dialog(tab="Advanced"));

      replaceable function pLoss =
          Chemical.Processes.Internal.FlowResistance.pleaseSelectChemicalPotentialLoss
        constrainedby
        Chemical.Processes.Internal.FlowResistance.partialChemicalPotentialLoss "ChemicalPotential loss function"
        annotation(choicesAllMatching=true, Documentation(info="<html>
<u>ChemicalPotential loss function used in the flow resistance.</u>
</html>"));

    protected
      Modelica.Units.SI.Density rho_rear_in=max(rho_min, Medium.density(rear.state_forwards)) "density of medium entering";
      Modelica.Units.SI.DynamicViscosity mu_rear_in=Medium.dynamicViscosity(rear.state_forwards) "dynamic viscosity of medium entering";

      Modelica.Units.SI.Density rho_fore_in=max(rho_min, Medium.density(fore.state_rearwards)) "density of medium entering";
      Modelica.Units.SI.DynamicViscosity mu_fore_in=Medium.dynamicViscosity(fore.state_rearwards) "dynamic viscosity of medium entering";

    equation
      //forwards model
      du_fore = -pLoss(n_flow, rho_rear_in, mu_rear_in, r, l);
      h_fore_out = h_rear_in;
      Xi_fore_out = Xi_rear_in;

      //rearwards model
      du_rear = -pLoss(-n_flow, rho_fore_in, mu_fore_in, r, l);
      h_rear_out = h_fore_in;
      Xi_rear_out = Xi_fore_in;

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Ellipse(
              extent={{-56,54},{64,-66}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{-100,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{40,0},{-48,0}},
              color={158,66,200},
              thickness=0.5,
              pattern=LinePattern.Dash),
            Line(
              points={{-44,-40},{0,-10},{44,-40}},
              color={158,66,200},
              thickness=0.5,
              smooth=Smooth.Bezier),
            Line(
              points={{-44,-15},{0,15},{44,-15}},
              color={158,66,200},
              thickness=0.5,
              smooth=Smooth.Bezier,
              origin={0,25},
              rotation=180)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Undirected implementation of the FlowResistance with different selectable flow resistance functions (laminar, laminar-turbulent, linear-quadratic). The output potential can be clipped to a certain value.</u>
</html>"));
    end FlowResistance;

    model TransportDelay "Delay chemical state depending on fluid speed"
      extends Chemical.Undirected.Interfaces.SISOBiFlow(final cliu_u_out=
            false);

      parameter Modelica.Units.SI.Length l "Length of Delay Pipe";
      parameter Modelica.Units.SI.Radius r "Radius of Delay Pipe";
      parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimal Density" annotation (Dialog(tab="Advanced"));
      parameter Real v_min(unit="1/s", min=0) = 0.01 "Minimum nondimensional speed"
        annotation(Dialog(tab="Advanced"));
      parameter Real v_max(unit="1/s", min=0) = 50 "Maximum nondimensional speed"
        annotation(Dialog(tab="Advanced"));

      constant Medium.ThermodynamicState state_0 = Chemical.Interfaces.SubstanceState(u=Medium.u_default,h= Medium.h_default);

      constant Modelica.Units.SI.SpecificVolume v_0=1/Medium.density(state_0);

      Real x(unit="1");
      Real v(unit="1/s");

    protected
      Modelica.Units.SI.SpecificInternalEnergy u_rear_in=Medium.specificInternalEnergy(rear.state_forwards);
      Modelica.Units.SI.SpecificInternalEnergy u_fore_in=Medium.specificInternalEnergy(fore.state_rearwards);
      Modelica.Units.SI.SpecificVolume v_rear_in=1/max(rho_min, Medium.density(rear.state_forwards));
      Modelica.Units.SI.SpecificVolume v_fore_in=1/max(rho_min, Medium.density(fore.state_rearwards));

      Modelica.Units.SI.SpecificInternalEnergy u_rear_out;
      Modelica.Units.SI.SpecificInternalEnergy u_fore_out;
      Modelica.Units.SI.SpecificVolume v_rear_out;
      Modelica.Units.SI.SpecificVolume v_fore_out;

      Modelica.Units.SI.Area A=r^2*Modelica.Constants.pi;

    initial equation
      x = 0;

    equation
      if n_flow >= 0 then
        v = min(v_max, max(v_min, n_flow*v_rear_in/A/l));
      else
        v = -min(v_max, max(v_min, -n_flow*v_fore_in/A/l));
      end if;

      der(x) = v;

      (u_rear_out,u_fore_out) = spatialDistribution(u_rear_in, u_fore_in,
        x, v>=0,
        initialPoints = {0.0,1.0},
        initialValues = {Medium.specificInternalEnergy(state_0), Medium.specificInternalEnergy(state_0)});
      (v_rear_out,v_fore_out) = spatialDistribution(v_rear_in, v_fore_in,
        x, v>=0,
        initialPoints = {0.0,1.0},
        initialValues = {v_0, v_0});

      for i in 1:Medium.nXi loop
        (Xi_rear_out[i], Xi_fore_out[i]) = spatialDistribution(Xi_rear_in[i], Xi_fore_in[i],
          x, v>=0,
          initialPoints = {0.0,0.0},
          initialValues = {Xi_0[i], Xi_0[i]});
      end for;

        //forwards model
      du_fore = 0;
      h_fore_out = u_fore_out + u_fore_out * v_fore_out;
      Xi_fore_out = Xi_rear_in;

      //rearwards model
      du_rear = 0;
      h_rear_out = u_rear_out + u_rear_out * v_rear_out;
      Xi_rear_out = Xi_fore_in;
      annotation (Documentation(info="<html>
<u>Undirected implementation of the transport delay.</u>
<u>Delays the temperature and massFraction, not potential, since potential differences propagate with speed of sound, and since delaying only steady state potential u not inertial potential r might lead to undesirable behavior.</u>
<u>Note that this component uses the spatialDistribution operator, that has some artefacts (see Fig. 1) for high and low non-dimensional speeds v (possibly due to inerpolation or extrapolation of the function). Therefore minimum and maximum speed in the non-dimensional coordinate x (inlet @ x=0, outlet @ x=1) is limited. The default limits are [0.01, 50], so the delay is limited by default to [0.02s, 100s]. This limit can be adjusted in the advanced parameters tab.</u>
<u><img src=\"modelica://Chemical/Resources/Doku/Chemical.Processes.Tests.TransportDelay_artefacts2.PNG\"/> <img src=\"modelica://Chemical/Resources/Doku/Chemical.Processes.Tests.TransportDelay_artefacts.PNG\"/> </u>
<u style=\"margin-left: 250px;\">Fig. 1: artefacts of the TransportDelay</u>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Ellipse(
              extent={{-56,54},{64,-66}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{-100,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
               points={{-25,-43},{25,-43},{-25,43},{25,43}},
               lineColor= {158,66,200},
               lineThickness=0.5),
            Line(
              points={{0,0},{0,-30}},
              color={158,66,200},
              thickness=0.5,
              pattern=LinePattern.Dot),
            Polygon(
               points={{-25,-43},{25,-43}, {0,-30}},
               lineColor= {158,66,200},
               fillColor={158,66,200},
               fillPattern=FillPattern.Solid),
            Polygon(
               points={{-15,26},{15,26}, {0,0}},
               lineColor= {158,66,200},
               fillColor={158,66,200},
               fillPattern=FillPattern.Solid)}));
    end TransportDelay;

    model ConductionElement "Volume with quasi-sationary mass and heatport"
      extends Internal.PartialConductionElement;

      parameter Boolean resistanceFromAU = true
        "= true, if thermal conductance given by U*A"
        annotation(Dialog(group="Thermal Conductance"));
      parameter Modelica.Units.SI.Area A=1 "Contact area of element with medium" annotation (Dialog(group="Thermal Conductance", enable=resistanceFromAU));
      parameter Modelica.Units.SI.CoefficientOfHeatTransfer U=200 "Heat transfer coefficient to medium"
        annotation (Dialog(group="Thermal Conductance", enable=resistanceFromAU));
      parameter Modelica.Units.SI.ThermalConductance k_par=200 "Thermal conductance heatport->fluid"
        annotation (Dialog(group="Thermal Conductance", enable=not resistanceFromAU));

    equation
      k = (if noEvent(resistanceFromAU) then U*A else k_par);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>
Undirected implementation of the Conduction Element.
</u>
<u>
This model is an element with a fixed volume (fig. 1). The mass in the volume is
assumed quasi-stationary (statically computed with volume and density), and the
fore massflow is coupled to the rear massflow.
</u>
<u>
<strong>Because of this the ConductionElement cannot be used as a loop breaker.</strong>
</u>
<u>
The advantage is that multiple ConductionElements can be put behind each other
without worrying about oscillations or fast eigenvalues between their masses.
The ConductionElement implements equations for conservation of mass and energy
for the fluid mass contained within it.
</u>
<u>
For further documentation see the documentation of the
<a href=\"modelica://Chemical.Undirected.Processes.Internal.PartialConductionElement\">motherclass</a>.
</u>
</html>"));
    end ConductionElement;

    package Tests "Tests for top level components of undirected"
      extends Modelica.Icons.ExamplesPackage;

      model TestFlowResistance "Test for the undirected flow resistance"
        extends Modelica.Icons.Example;

        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(substanceData=Chemical.Substances.Water_liquid(),
          u0_par=100000) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-30,0})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(
          potentialFromInput=true,
          u0_par=110000) annotation (Placement(transformation(extent={{20,-10},{40,10}})));
        inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=0.01) annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
        Modelica.Blocks.Sources.Step step(
          height=-80000,
          offset=140000,
          startTime=5)
          annotation (Placement(transformation(extent={{60,-6},{48,6}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear1(substanceData=Chemical.Substances.Water_liquid(),
          potentialFromInput=true)
                         annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-28,-38})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore1(
          potentialFromInput=false,
          u0_par=100000) annotation (Placement(transformation(extent={{22,-48},{42,-28}})));
        Reaction reaction(
          productsSubstanceData={Chemical.Substances.Water_liquid()},
          nS=1,
          nP=1) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Reaction reaction1(
          productsSubstanceData={Chemical.Substances.Water_liquid()},
          nS=1,
          nP=1) annotation (Placement(transformation(extent={{-8,-48},{12,-28}})));
      equation
        connect(step.y, boundary_fore.u0_var)
          annotation (Line(points={{47.4,0},{40,0},{40,6},{32,6}},
                                                         color={0,0,127}));
        connect(boundary_rear1.u0_var, boundary_fore.u0_var)
          annotation (Line(points={{-30,-44},{-40,-44},{-40,-20},{40,-20},{40,6},{32,6}}, color={0,0,127}));
        connect(boundary_rear.fore, reaction.substrates[1]) annotation (Line(
            points={{-20,0},{-10,0}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction.products[1], boundary_fore.rear) annotation (Line(
            points={{10,0},{20,0}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear1.fore, reaction1.substrates[1]) annotation (Line(
            points={{-18,-38},{-8,-38}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction1.products[1], boundary_fore1.rear) annotation (Line(
            points={{12,-38},{22,-38}},
            color={158,66,200},
            thickness=0.5));
        annotation (
          Icon(graphics,
               coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
          Documentation(info="<html>
<u>Test for the undirected flow resistance.</u>
<u><br>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
      end TestFlowResistance;

      model SimpleReaction "The simple chemical reaction A<->B with equilibrium B/A = 2"
         extends Modelica.Icons.Example;

        constant Real K = 2 "Dissociation constant of the reaction";

        constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
        constant Real R = Modelica.Constants.R "Gas constant";

        Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

        Boundaries.Substance          A(
          useRear=false,
          useFore=true,
          useSolution=true,             use_mass_start=false, amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-52,-8},{-32,12}})));

        Reaction reaction2_1(
          productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(K))},
          nS=1,
          nP=1) annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
        Boundaries.Substance B(
          useRear=true,
          useSolution=true,
          use_mass_start=false,
          amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{42,-8},{62,12}})));

        inner Modelica.Fluid.System system annotation (Placement(transformation(extent={{58,64},{78,84}})));
      equation
        connect(A.solution, solution.solution) annotation (Line(
            points={{-48,-8},{-48,-92},{60,-92},{60,-98}},
            color={127,127,0}));
        connect(B.solution, solution.solution) annotation (Line(points={{46,-8},{46,-92},{60,-92},{60,-98}},
                                           color={127,127,0}));
        connect(A.fore, reaction2_1.substrates[1]) annotation (Line(
            points={{-32,2},{-10,2}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction2_1.products[1], B.rear) annotation (Line(
            points={{10,2},{42,2}},
            color={158,66,200},
            thickness=0.5));
        annotation (Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Simple reaction demonstrating equilibria between substance A and substance B, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that mole fraction (A.x and B.x) are always summed to 1 for the solution.</p>
</html>"),experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
      end SimpleReaction;

      model SimpleReaction2 "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
         extends Modelica.Icons.Example;

        constant Real Kb(unit="kg/mol") = 2
          "Molarity based dissociation constant of the reaction with one more reactant";

        constant Real Kx(unit="1") = Kb*55.508
          "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

        constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
        constant Real R = Modelica.Constants.R "Gas constant";

        Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

        Boundaries.Substance          A(
          useFore=true,
          useSolution=true,
          substanceData(MolarWeight=1),
          use_mass_start=false,
          amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,4},{-14,24}})));
        Reaction reaction2_1(
          productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(MolarWeight=2, DfG=-R*T_25degC*log(Kx))},
          nS=2,
          nP=1) annotation (Placement(transformation(extent={{4,-8},{24,12}})));
        Boundaries.Substance          B(
          useFore=true,
          useSolution=true,
          substanceData(MolarWeight=1),
          use_mass_start=false,
          amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
        Boundaries.Substance C(
          useRear=true,
          useSolution=true,
          use_mass_start=false,
          amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{48,-8},{68,12}})));

      equation
        connect(A.solution, solution.solution) annotation (Line(
            points={{-30,4},{-30,-90},{60,-90},{60,-98}},
            color={127,127,0}));
        connect(C.solution, solution.solution) annotation (Line(points={{52,-8},{66,-8},{66,-90},{60,-90},{60,-98}},
                                                   color={127,127,0}));
        connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},
              {-30,-90},{60,-90},{60,-98}},color={127,127,0}));

        connect(A.fore, reaction2_1.substrates[1])
          annotation (Line(
            points={{-14,14},{-4,14},{-4,1.75},{4,1.75}},
            color={158,66,200},
            thickness=0.5));
        connect(B.fore, reaction2_1.substrates[2])
          annotation (Line(
            points={{-14,-14},{-8,-14},{-8,2.25},{4,2.25}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction2_1.products[1], C.rear) annotation (Line(
            points={{24,2},{48,2}},
            color={158,66,200},
            thickness=0.5));
        annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
      end SimpleReaction2;

      model SimpleReaction22 "The simple chemical reaction A+B<->C+D with equilibrium [C]*[D]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
         extends Modelica.Icons.Example;

        constant Real Kb(unit="kg/mol") = 2
          "Molarity based dissociation constant of the reaction with one more reactant";

        constant Real Kx(unit="1") = Kb*55.508
          "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

        constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
        constant Real R = Modelica.Constants.R "Gas constant";

        Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

        Chemical.Boundaries.Substance A(use_mass_start=false, amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
        Chemical.Processes.Reaction reaction2_1(
          productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(Kx)),
              Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(Kx))},
          nS=2,
          nP=2) annotation (Placement(transformation(extent={{4,-8},{24,12}})));
        Chemical.Boundaries.Substance B(use_mass_start=false, amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
        Chemical.Boundaries.Substance C(amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{48,-8},{68,12}})));

        Chemical.Boundaries.Substance D(amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{44,-34},{64,-14}})));
        inner DropOfCommons dropOfCommons(L=1e-3) annotation (Placement(transformation(extent={{52,56},{72,76}})));
      equation
        connect(A.solution, solution.solution) annotation (Line(
            points={{-30,2},{-30,-90},{60,-90},{60,-98}},
            color={127,127,0}));
        connect(C.solution, solution.solution) annotation (Line(points={{52,-8},{66,-8},{66,-90},{60,-90},{60,-98}},
                                                   color={127,127,0}));
        connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},
              {-30,-90},{60,-90},{60,-98}},color={127,127,0}));

        connect(B.outlet, reaction2_1.substrates[1]) annotation (Line(
            points={{-14,-14},{-4,-14},{-4,1.75},{4,1.75}},
            color={158,66,200},
            thickness=0.5));
        connect(A.outlet, reaction2_1.substrates[2]) annotation (Line(
            points={{-14,12},{-4,12},{-4,2.25},{4,2.25}},
            color={158,66,200},
            thickness=0.5));
        connect(D.solution, solution.solution) annotation (Line(points={{48,-34},{60,-34},{60,-98}}, color={127,127,0}));
        connect(reaction2_1.products[1], C.inlet) annotation (Line(
            points={{24,1.75},{36,1.75},{36,2},{48,2}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction2_1.products[2], D.inlet) annotation (Line(
            points={{24,2.25},{34,2.25},{34,-24},{44,-24}},
            color={158,66,200},
            thickness=0.5));
        annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),experiment(StopTime=100, __Dymola_Algorithm="Dassl"));
      end SimpleReaction22;

      model ExothermicReaction "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"

         extends Modelica.Icons.Example;

        parameter Modelica.Units.SI.MolarEnergy ReactionEnthalpy=-55000;

        Chemical.Solution thermal_isolated_solution(useMechanicPorts=true, ConstantTemperature=false)
          annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));
        Boundaries.Substance          A(
          useFore=true,
          useSolution=true,             use_mass_start=false, amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
        Reaction reaction2_2(
          productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfH=ReactionEnthalpy)},
          nS=1,
          nP=1) annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
        Boundaries.Substance B(
          useRear=true,
          useSolution=true,
          use_mass_start=false,
          amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{20,-60},{40,-40}})));

        Chemical.Solution solution_at_constant_temperature(useMechanicPorts=true, useThermalPort=true)
          annotation (Placement(transformation(extent={{-100,0},{98,94}})));
        Boundaries.Substance          A1(
          useFore=true,
          useSolution=true,              use_mass_start=false, amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
        Reaction reaction2_1(
          productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfH=ReactionEnthalpy)},
          nS=1,
          nP=1) annotation (Placement(transformation(extent={{-8,40},{12,60}})));
        Boundaries.Substance B1(
          useRear=true,
          useSolution=true,
          use_mass_start=false,
          amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{20,40},{40,60}})));

        //  Modelica.SIunits.HeatFlowRate q
        //    "Heat flow to environment to reach constant temperature";
        Modelica.Units.SI.Temperature t
          "Temperature if the solution is ideally thermal isolated from environment";
        Chemical.Boundaries.Substance H2O(substanceData=Chemical.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{20,4},{40,24}})));
        Chemical.Boundaries.Substance H2O1(substanceData=Chemical.Substances.Water_liquid(), mass_start=1)
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
            points={{24,-60},{24,-64},{58.4,-64},{58.4,-99.06}},
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
            points={{12,-50},{20,-50}},
            color={158,66,200},
            thickness=0.5));
        annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Demonstration of exotermic reaction with perfect cooling (i.e. connected fixed temperature to the HeatPort) and thermally insulated (HetPort unconnected). See solution_(...).T</p>
</html>"),experiment(StopTime=10, __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}})));
      end ExothermicReaction;

      model ConductionElement "Test for ConductionElement"
        extends Modelica.Icons.Example;

        replaceable package Medium =
            Chemical.Media.myMedia.Incompressible.Examples.Glycol47                          constrainedby
          Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium Model"
          annotation(choicesAllMatching=true, Documentation(info="<html>
<u>
Medium Model for the test. Be aware that the Component is mainly
meant for liquids with low compressablility.
</u>
</html>"));

        Chemical.Undirected.Processes.ConductionElement conductionElement(
          redeclare package stateOfMatter = stateOfMatter,
          L=0.2,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=-1,
          V(displayUnit="l") = 0.001,
          A=35,
          U=500,
          init=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.T,
          T_0=263.15) annotation (Placement(transformation(extent={{-10,60},{10,80}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(
          redeclare package stateOfMatter = stateOfMatter,
          T0_par=328.15,
          u0_par=100000) annotation (Placement(transformation(extent={{20,60},{40,80}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par=288.15) annotation (Placement(transformation(extent={{-40,60},{-20,80}})));
        Modelica.Blocks.Sources.Step step(
          height=2e2,
          offset=0.999e5,
          startTime=0.33)
          annotation (Placement(transformation(extent={{-70,60},{-50,80}})));
        Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=373.15)
          annotation (Placement(transformation(extent={{80,50},{60,70}})));
        inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{60,-70},{80,-50}})));
        Chemical.Undirected.Processes.ConductionElement conductionElement1(
          redeclare package stateOfMatter = stateOfMatter,
          L=0.2,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          V(displayUnit="l") = 0.001,
          A=35,
          U=500,
          init=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.h,
          T_0=263.15,
          h_0=1000,
          neglectChemicalPotentialChanges=false) annotation (Placement(transformation(extent={{-10,30},{10,50}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore1(
          redeclare package stateOfMatter = stateOfMatter,
          T0_par=328.15,
          u0_par=100000) annotation (Placement(transformation(extent={{20,30},{40,50}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear1(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par=288.15) annotation (Placement(transformation(extent={{-40,30},{-20,50}})));
        Modelica.Blocks.Sources.Ramp ramp(
          height=2e2,
          duration=0.001,
          offset=0.999e5,
          startTime=0.33)
          annotation (Placement(transformation(extent={{-70,30},{-50,50}})));
        Chemical.Undirected.Processes.ConductionElement conductionElement2(
          redeclare package stateOfMatter = stateOfMatter,
          L=0.2,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=-1,
          V(displayUnit="l") = 0.001,
          enforce_global_energy_conservation=true,
          A=35,
          U=500,
          init=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.fore,
          T_0=263.15) annotation (Placement(transformation(extent={{-10,0},{10,20}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore2(
          redeclare package stateOfMatter = stateOfMatter,
          T0_par=328.15,
          u0_par=100000) annotation (Placement(transformation(extent={{20,0},{40,20}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear2(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par=288.15) annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
        Modelica.Blocks.Sources.Step step1(
          height=2e2,
          offset=0.999e5,
          startTime=0.33)
          annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
        Chemical.Undirected.Processes.ConductionElement conductionElement3(
          redeclare package stateOfMatter = stateOfMatter,
          L=0.2,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=-1,
          V(displayUnit="l") = 0.001,
          A=35,
          U=500,
          init=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.rear,
          T_0=263.15) annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore3(
          redeclare package stateOfMatter = stateOfMatter,
          T0_par=328.15,
          u0_par=100000) annotation (Placement(transformation(extent={{20,-30},{40,-10}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear3(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par=288.15) annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
        Modelica.Blocks.Sources.Step step2(
          height=2e2,
          offset=0.999e5,
          startTime=0.33)
          annotation (Placement(transformation(extent={{-70,-30},{-50,-10}})));
        Chemical.Undirected.Processes.ConductionElement conductionElement4(
          redeclare package stateOfMatter = stateOfMatter,
          L=0.2,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=0,
          V(displayUnit="l") = 0.001,
          enforce_global_energy_conservation=true,
          A=35,
          U=500,
          init=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.port,
          T_0=263.15) annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore4(
          redeclare package stateOfMatter = stateOfMatter,
          T0_par=328.15,
          u0_par=100000) annotation (Placement(transformation(extent={{20,-60},{40,-40}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear4(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par=288.15) annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
        Modelica.Blocks.Sources.Step step3(
          height=2e2,
          offset=0.999e5,
          startTime=0.33)
          annotation (Placement(transformation(extent={{-70,-60},{-50,-40}})));
        Chemical.Undirected.Processes.ConductionElement conductionElement5(
          redeclare package stateOfMatter = stateOfMatter,
          L=0.2,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=-1,
          V(displayUnit="l") = 0.001,
          A=35,
          U=500,
          init=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.port,
          T_0=263.15) annotation (Placement(transformation(extent={{-10,-90},{10,-70}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore5(
          redeclare package stateOfMatter = stateOfMatter,
          T0_par=328.15,
          u0_par=100000) annotation (Placement(transformation(extent={{20,-90},{40,-70}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear5(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par=288.15) annotation (Placement(transformation(extent={{-40,-90},{-20,-70}})));
        Modelica.Blocks.Sources.Step step4(
          height=2e2,
          offset=0.999e5,
          startTime=0.33)
          annotation (Placement(transformation(extent={{-70,-90},{-50,-70}})));
      equation
        connect(step.y, boundary_rear.u0_var) annotation (Line(points={{-49,70},{-40,70},{-40,76},{-32,76}},
                                                color={0,0,127}));
        connect(conductionElement.rear, boundary_rear.fore) annotation (Line(
            points={{-10,70},{-20,70}},
            color={158,66,200},
            thickness=0.5));
        connect(conductionElement.fore, boundary_fore.rear) annotation (Line(
            points={{10,70},{20,70}},
            color={158,66,200},
            thickness=0.5));
        connect(ramp.y, boundary_rear1.u0_var) annotation (Line(points={{-49,40},{-40,40},{-40,46},{-32,46}},
                                                                                              color={0,0,127}));
        connect(conductionElement1.rear, boundary_rear1.fore)
          annotation (Line(
            points={{-10,40},{-20,40}},
            color={158,66,200},
            thickness=0.5));
        connect(conductionElement1.fore, boundary_fore1.rear) annotation (Line(
            points={{10,40},{20,40}},
            color={158,66,200},
            thickness=0.5));
        connect(step1.y, boundary_rear2.u0_var) annotation (Line(points={{-49,10},{-40,10},{-40,16},{-32,16}},
                                                                                             color={0,0,127}));
        connect(conductionElement2.rear, boundary_rear2.fore) annotation (Line(
            points={{-10,10},{-20,10}},
            color={158,66,200},
            thickness=0.5));
        connect(conductionElement2.fore, boundary_fore2.rear) annotation (Line(
            points={{10,10},{20,10}},
            color={158,66,200},
            thickness=0.5));
        connect(step2.y, boundary_rear3.u0_var) annotation (Line(points={{-49,-20},{-40,-20},{-40,-14},{-32,-14}},
                                                                                               color={0,0,127}));
        connect(conductionElement3.rear, boundary_rear3.fore)
          annotation (Line(
            points={{-10,-20},{-20,-20}},
            color={158,66,200},
            thickness=0.5));
        connect(conductionElement3.fore, boundary_fore3.rear) annotation (Line(
            points={{10,-20},{20,-20}},
            color={158,66,200},
            thickness=0.5));
        connect(step3.y, boundary_rear4.u0_var) annotation (Line(points={{-49,-50},{-40,-50},{-40,-44},{-32,-44}},
                                                                                               color={0,0,127}));
        connect(conductionElement4.rear, boundary_rear4.fore)
          annotation (Line(
            points={{-10,-50},{-20,-50}},
            color={158,66,200},
            thickness=0.5));
        connect(conductionElement4.fore, boundary_fore4.rear) annotation (Line(
            points={{10,-50},{20,-50}},
            color={158,66,200},
            thickness=0.5));
        connect(step4.y, boundary_rear5.u0_var) annotation (Line(points={{-49,-80},{-40,-80},{-40,-74},{-32,-74}},
                                                                                               color={0,0,127}));
        connect(conductionElement5.rear, boundary_rear5.fore)
          annotation (Line(
            points={{-10,-80},{-20,-80}},
            color={158,66,200},
            thickness=0.5));
        connect(conductionElement5.fore, boundary_fore5.rear) annotation (Line(
            points={{10,-80},{20,-80}},
            color={158,66,200},
            thickness=0.5));
        connect(fixedTemperature.port, conductionElement1.heatPort) annotation (Line(points={{60,60},{40,60},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
        connect(conductionElement.heatPort, conductionElement1.heatPort)
          annotation (Line(points={{0,79.8},{0,88},{40,88},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
        connect(conductionElement2.heatPort, conductionElement1.heatPort)
          annotation (Line(points={{0,19.8},{0,26},{40,26},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
        connect(conductionElement3.heatPort, conductionElement1.heatPort)
          annotation (Line(points={{0,-10.2},{0,-4},{40,-4},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
        connect(conductionElement4.heatPort, conductionElement1.heatPort)
          annotation (Line(points={{0,-40.2},{0,-34},{40,-34},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
        connect(conductionElement5.heatPort, conductionElement1.heatPort)
          annotation (Line(points={{0,-70.2},{0,-66},{40,-66},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
        annotation (experiment(StopTime=1, Tolerance=1e-6, Interval=0.001),
        Documentation(info="<html>
<u>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
      end ConductionElement;

      model TransportDelay "Test for transport delay"
        extends Modelica.Icons.Example;

        replaceable package Medium = Chemical.Media.myMedia.Air.DryAirNasa constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
                                                               "Medium Model"
          annotation (Documentation(info="<html>
<u>
Medium model for the test. Can be anything.
</u>
</html>"));

        Chemical.Processes.TransportDelay transportDelay(
          redeclare package stateOfMatter = stateOfMatter,
          l=100,
          r(displayUnit="mm") = 0.015)
          annotation (Placement(transformation(extent={{30,30},{50,50}})));
        Chemical.Processes.FlowResistance flowResistance(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=100,
          l(displayUnit="mm") = 0.008,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (
            k=1e4))
          annotation (Placement(transformation(extent={{-10,30},{10,50}})));

        Chemical.Boundaries.Sink sink(redeclare package stateOfMatter =
              stateOfMatter,
            u0_par=100000)
          annotation (Placement(transformation(extent={{70,30},{90,50}})));
        Chemical.Boundaries.Source source(
          redeclare package stateOfMatter = stateOfMatter,
          temperatureFromInput=true,
          potentialFromInput=true)
          annotation (Placement(transformation(extent={{-48,30},{-28,50}})));
        Modelica.Blocks.Sources.Trapezoid
                                     trapezoid(
          amplitude=-6e3,
          rising=0.2,
          width=0.2,
          falling=0.2,
          period=0.8,
          offset=1.04e5,
          startTime=0.2)
          annotation (Placement(transformation(extent={{-94,10},{-74,30}})));
        Modelica.Blocks.Sources.Ramp ramp1(
          height=20,
          duration=1.5,
          offset=273,
          startTime=0.3)
          annotation (Placement(transformation(extent={{-94,50},{-74,70}})));
        inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{60,-80},{80,-60}})));

        Chemical.Undirected.Processes.TransportDelay transportDelay1(
          redeclare package stateOfMatter = stateOfMatter,
          l=100,
          r(displayUnit="mm") = 0.015) annotation (Placement(transformation(extent={{30,-50},{50,-30}})));
        FlowResistance flowResistance1(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=100,
          l(displayUnit="mm") = 0.008,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (
            k=1e4))
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=100000)
          annotation (Placement(transformation(extent={{70,-50},{90,-30}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package stateOfMatter = stateOfMatter,
          temperatureFromInput=true,
          potentialFromInput=true) annotation (Placement(transformation(extent={{-48,-50},{-28,-30}})));
        Modelica.Blocks.Sources.Trapezoid trapezoid1(
          amplitude=-6e3,
          rising=0.2,
          width=0.2,
          falling=0.2,
          period=0.8,
          offset=1.04e5,
          startTime=0.2)
          annotation (Placement(transformation(extent={{-94,-70},{-74,-50}})));
        Modelica.Blocks.Sources.Ramp ramp2(
          height=20,
          duration=1.5,
          offset=273,
          startTime=0.3)
          annotation (Placement(transformation(extent={{-94,-30},{-74,-10}})));
      equation
        connect(ramp1.y, source.T0_var) annotation (Line(points={{-73,60},{-60,60},{-60,40},{-40,40}},
                                 color={0,0,127}));
        connect(trapezoid.y, source.u0_var) annotation (Line(points={{-73,20},{-60,20},{-60,46},{-40,46}},
                                        color={0,0,127}));
        connect(flowResistance.inlet, source.outlet) annotation (Line(
            points={{-10,40},{-28,40}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance.outlet, transportDelay.inlet) annotation (Line(
            points={{10,40},{30,40}},
            color={158,66,200},
            thickness=0.5));
        connect(transportDelay.outlet, sink.inlet) annotation (Line(
            points={{50,40},{70,40}},
            color={158,66,200},
            thickness=0.5));
        connect(ramp2.y, boundary_rear.T0_var) annotation (Line(points={{-73,-20},{-60,-20},{-60,-40},{-40,-40}},
                                                color={0,0,127}));
        connect(trapezoid1.y, boundary_rear.u0_var) annotation (Line(points={{-73,-60},{-60,-60},{-60,-34},{-40,-34}},
                                                     color={0,0,127}));
        connect(flowResistance1.rear, boundary_rear.fore) annotation (Line(
            points={{-10,-40},{-28,-40}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance1.fore, transportDelay1.rear) annotation (Line(
            points={{10,-40},{30,-40}},
            color={158,66,200},
            thickness=0.5));
        connect(transportDelay1.fore, boundary_fore.rear) annotation (Line(
            points={{50,-40},{62,-40},{62,-40},{70,-40}},
            color={158,66,200},
            thickness=0.5));
        annotation (
          experiment(
            StopTime=2.5,
            Tolerance=1e-6,
            Interval=0.0025,
            __Dymola_Algorithm="Dassl"),
          Icon(coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<u>Test for use Transport Delay to delay the thermodynamic state of the flow. </u>
<u>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
      end TransportDelay;

      model EnzymeKinetics "Basic enzyme kinetics"
        extends Modelica.Icons.Example;

        Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

        //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
        Boundaries.Substance P(
          useRear=true,
          useSolution=true,
          use_mass_start=false,amountOfSubstance_start=1e-8) annotation (Placement(transformation(extent={{72,-12},{92,8}})));

        Boundaries.Substance          S(
          useFore=true,
          useSolution=true,             use_mass_start=false, amountOfSubstance_start=100) annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));

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

        Boundaries.Substance ES(
          useRear=true,
          useFore=true,
          useSolution=true,
          initAmount=Chemical.Utilities.Types.InitializationMethods.state,
          use_mass_start=false, amountOfSubstance_start=tE/2) annotation (Placement(transformation(extent={{-8,-10},{12,10}})));
        Boundaries.Substance E(
          useRear=true,
          useFore=true,
          useSolution=true,
          use_mass_start=false,amountOfSubstance_start=tE/2) annotation (Placement(transformation(extent={{12,36},{-8,56}})));
        Reaction                    chemicalReaction(
          k_forward=1,
          productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-Modelica.Constants.R*298.15*log(2/Km))},
                nS=2,
          nP=1) annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));

        Processes.ForwardReaction chemicalReaction1(
          k_forward=k_cat,
          productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-Modelica.Constants.R*298.15*50),
              Chemical.Interfaces.Incompressible.SubstanceDataParameters()},
          nS=1,
          nP=2) annotation (Placement(transformation(extent={{24,-8},{44,12}})));

        Boundaries.Substance liquidWater(
          useSolution=true,              substanceData=Chemical.Substances.Water_liquid(),
          use_mass_start=true,                                                             mass_start=1)
          annotation (Placement(transformation(extent={{42,-80},{62,-60}})));
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
            points={{-22,0},{-8,0}},
            color={158,66,200},
            thickness=0.5));
        connect(ES.fore, chemicalReaction1.substrates[1])
          annotation (Line(
            points={{12,0},{18,0},{18,2},{24,2}},
            color={158,66,200},
            thickness=0.5));
        connect(S.fore, chemicalReaction.substrates[1])
          annotation (Line(
            points={{-72,-4},{-52,-4},{-52,-2},{-42,-2},{-42,-0.25}},
            color={158,66,200},
            thickness=0.5));
        connect(E.fore, chemicalReaction.substrates[2])
          annotation (Line(
            points={{-8,46},{-52,46},{-52,0.25},{-42,0.25}},
            color={158,66,200},
            thickness=0.5));
        connect(chemicalReaction1.products[1], P.rear)
          annotation (Line(
            points={{44,1.75},{58,1.75},{58,-2},{72,-2}},
            color={158,66,200},
            thickness=0.5));
        connect(E.rear, chemicalReaction1.products[2])
          annotation (Line(
            points={{12,46},{48,46},{48,48},{56,48},{56,2.25},{44,2.25}},
            color={158,66,200},
            thickness=0.5));
            annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",       info="<html>
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
</html>"),experiment(StopTime=100000, __Dymola_Algorithm="Dassl"),
          __Dymola_experimentSetupOutput);
      end EnzymeKinetics;
      annotation (Documentation(info="<html>
<u>Tests for top level components of the undirected chemical simulation package.</u>
</html>"));
    end Tests;

    package Internal "Internals package for Processes"
      extends Modelica.Icons.InternalPackage;

      partial model PartialConductionElement "Partial volume with quasisationary mass and heatport and undetermined heat transfer coefficient"
        extends Chemical.Undirected.Interfaces.SISOBiFlow(
                                      final cliu_u_out=false);

        parameter Modelica.Units.SI.Volume V=1 "Volume of the element";
        parameter Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement init=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.rear
          "Initialization method for h" annotation (Dialog(tab="Initialization", group="Enthalpy"));
        parameter Modelica.Units.SI.Temperature T_0=Medium.T_default "Initial Temperature" annotation (Dialog(
            tab="Initialization",
            group="Enthalpy",
            enable=(init == Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.T)));
        parameter Modelica.Units.SI.MolarEnthalpy h_0=Medium.h_default "Initial molar enthalpy" annotation (Dialog(
            tab="Initialization",
            group="Enthalpy",
            enable=(init == Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.h)));
        parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimal density" annotation (Dialog(tab="Advanced"));
        parameter Boolean neglectChemicalPotentialChanges = true "Neglect potential changes in energy equation"
          annotation(Dialog(tab="Advanced"));
        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg
          "Regularization massflow to switch between positive- and negative-massflow model" annotation (Dialog(tab="Advanced"));
        parameter Boolean enforce_global_energy_conservation = false "If true, exact global energy conservation is enforced by feeding back all energy stored locally back in the system"
          annotation(Dialog(tab="Advanced"));
        parameter Modelica.Units.SI.Time T_e=100 "Factor for feeding back energy."
          annotation (Dialog(tab="Advanced", enable=enforce_global_energy_conservation));

        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort(Q_flow=Q_flow, T=T_heatPort)
          annotation (Placement(transformation(extent={{-10,88},{10,108}})));

        Modelica.Units.SI.MolarEnthalpy h(start=Medium.h_default, stateSelect=StateSelect.prefer);

        Medium.ThermodynamicState state = Chemical.Interfaces.SubstanceState(u=u,h= h);
        Modelica.Units.SI.Temperature T=Medium.temperature(state);
        Modelica.Units.SI.ThermalConductance k "Thermal conductance heatport->fluid";

        Modelica.Units.SI.Energy deltaE_system(start=0, fixed=true) "Energy difference between n_flow*(h_in-h_out) and Q_flow";

        Modelica.Units.SI.Mass M;

      protected
        Modelica.Units.SI.ChemicalPotential u=Chemical.Undirected.Internal.regStep(
                  n_flow,
                  u_rear_in,
                  u_fore_in,
                  n_flow_reg);

        //h_in only in rhs of ODE--> h still smooth, better results at low massflow than using regStep
        Modelica.Units.SI.MolarEnthalpy h_in=if n_flow >= 0 then h_rear_in else h_fore_in;

        Modelica.Units.SI.Density rho=max(rho_min, Medium.density(state));
        Modelica.Units.SI.MolarEnthalpy h_out=state.h;

        Modelica.Units.SI.Temperature T_heatPort;
        Modelica.Units.SI.HeatFlowRate Q_flow;

      initial equation
        if init == Internal.InitializationMethodsCondElement.T then
           h = Medium.setState_pTX(u, T_0, Xi);
        elseif init == Internal.InitializationMethodsCondElement.h then
           h = h_0;
        elseif init == InitializationMethodsCondElement.port then
           h = h_in;
        elseif init == Internal.InitializationMethodsCondElement.fore then
           h = h_fore_in;
         else
           h = h_rear_in;
         end if;

      equation
        //mass is assumed to be quasisationary-> this violates the conservation of mass since n_flow_in = -n_flow_out. see documentation/information
        M = V*rho;

        // Mass is assumed to be (quasi-)stationary-> assumption: der(M)=0.
        // V=const
        // lhs: der(U) = der(H - pV) = der(M)*h +  M*der(h) - V*der(u) = M*der(h)
        // one can not simply add der(M)*h to the left side, since the right side is also assuming sationary mass (n_flow_in=-n_flow_out)
        // on RHS: Q_flow + n_flow_in*h_in + n_flow_out*h = Q_flow + n_flow_in*h_in + (n_flow_in-n_flow_in+ n_flow_out)*h = Q_flow + n_flow_in*(h_in-h) + der(M)*h = Q_flow + n_flow*(h_in-h)
        // or:     Q_flow + n_flow_in*h_in + n_flow_out*h = Q_flow + (n_flow_in+n_flow_out-n_flow_out)*h_in + n_flow_out*h = Q_flow + n_flow_out*(h-h_in) + der(M)*h_in = Q_flow + n_flow*(h_in - h)
        // so there is a resulting term on the RHS with der(M)*(h(k) + h_out(1-k)) with k in [0,1]; that also would have to be accounted for when der(M) not 0 on LHS
        // to improve the equation u_in and M would have to be filtered and a value for k would have to be found, all that resulted in only moderate improvements
        // therefore energy that is accumulated with system border around the component (deltaE_system) is slowly fed back into the system when enforce_global_energy_conservation is true
        //(then energy stored in the system e.g. by evaproation/condension, while still being visible short-term, neglegted in logterm)
        if not neglectChemicalPotentialChanges then
          M*der(h) = Q_flow + abs(n_flow)*(h_in - h_out) + V*der(u) + (if enforce_global_energy_conservation then deltaE_system/T_e else 0);
        else
          //neglegt V*der(u), since u might not be smooth -> to notice a difference der(u) must be about 1e7 Pa/s. see documentation/information
          M*der(h) = Q_flow + abs(n_flow)*(h_in - h_out) + (if enforce_global_energy_conservation then deltaE_system/T_e else 0);
        end if;
        if enforce_global_energy_conservation then
          der(deltaE_system) = Q_flow + abs(n_flow)*(h_in - h_out);
        else
          deltaE_system = 0;
        end if;

        Q_flow = k*(T_heatPort - T);

        //forwards model
        du_fore = 0;
        h_fore_out = h;
        Xi_fore_out = Xi_rear_in;

        //rearwards model
        du_rear = 0;
        h_rear_out = h;
        Xi_rear_out = Xi_fore_in;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
             Line(
               points={{-100,0},{100,0}},
               thickness=0.5,
               color={158,66,200}),
             Ellipse(
               extent={{-70,-70},{70,70}},
               lineColor={158,66,200},
               lineThickness=0.5,
               fillColor={194,138,221},
               fillPattern=FillPattern.Solid,
               pattern=LinePattern.Solid),
             Line(
               points={{-50,-30},{50,-30}},
               color={238,46,47}),
             Line(
               points={{-50,-15},{50,-15}},
               color={238,46,47}),
             Line(
               points={{-50,0},{50,0}},
               color={238,46,47}),
             Line(
               points={{-50,15},{50,15}},
               color={238,46,47}),
             Line(
               points={{-50,30},{50,30}},
               color={238,46,47}),
             Line(
               points={{0,100},{0,-30}},
               color={238,46,47})}),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<u>
Undirected implementation of the Conduction Element.
</u>
<u>
This model is an element with a fixed volume (fig. 1). The mass in the volume is 
assumed quasi-stationary (statically computed with volume and density), and the 
fore massflow is coupled to the rear massflow. <strong>Because of this the 
ConductionElement cannot be used as a loop breaker</strong>.
The advantage is that multiple ConductionElements can be put behind each other 
without worrying about oscillations or fast eigenvalues between their masses. 
The ConductionElement implements equations for conservation of mass and energy 
for the fluid mass contained within it.
</u>
<u>
Different to the unidirectional ConductionElement, the model for forward massflow
(see fig. 1) is valid for both flow directions.
</u>
<u>
Additionally the undirected Conduction Element offers more initialization methods 
for h. Additionally to initialize T or h by a parameter, one can choose to 
initialize with the incoming enthalpy from either one of the two ports, or use the 
correct one, depending on the massflow (option &apos;port&apos;). The last option 
can lead to large nonlinear initialization problems, we advice to choose the port 
to initialize h from if known in advance (options &apos;rear&apos; or &apos;fore&apos;).
</u>
<u>
The ConductionElement makes different assumptions:
</u>
<ul>
  <li>
    Quasistationary Mass:
    <br>
    n_flow_rear = - n_flow_fore &amp; M = rho * V
    (this assumption violates the conservation of mass for changing densities,
    since the mass in the element can change although inflow and outflow are the
    same)
    <br>
    der(H) = der(M*h) = M*der(h)
    (This assumption violates the conservation of energy for changing densities,
    since then the mass M of fluid in the element is no longer constant)
  </li>
  <li>
    Neglection of der(u) in the energy equation
    <br>
    V*der(u) = 0
    (this assumption violates the conservation of energy for changing potentials.
    For a noticeable difference in the testcase the der(u) must be in the
    order of 1e5 Pa/s).
    <br>
    This assumption can be turned off by setting neglectChemicalPotentialChanges=false
    (true by default) in the Advanced tab. <strong>This option requires the 
    fore and rear input potentials to be smooth.</strong>
  </li>
</ul>
<u>
Due to these assumptions minor violations in the global energy conservation can 
occur.
With the flag enforce_global_energy_conservation in the &quot;Advanced&quot; tab 
is set true (Default: false), long-term energy storage in the ConductionElement 
is sacrificed to hold global energy conservation.
</u>
<div>
<img src=\"modelica://Chemical/Resources/Doku/Chemical.Processes.ConductionElement_positive.png\"/>
</div>
<u>
fig. 1: positive massflow model
</u>
</html>"));
      end PartialConductionElement;

      type InitializationMethodsCondElement = enumeration(
          T "Temperature T_0",
          h "molar enthalpy h_0",
          fore "input state from fore",
          rear "input state from rear",
          port "regularized input state from fore or rear, depending on massflow") "Choices for initialization of state h of undirected ConductionElement"
        annotation (
          Icon(coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<u>
Choices for initialization of a state h.
</u>
</html>"));
      partial model PartialReactionWithSubstanceData "Chemical Reaction"
        extends Chemical.Undirected.Processes.Internal.PartialReaction;

        parameter stateOfMatter.SubstanceDataParameters productsSubstanceData[nP]
         annotation (choicesAllMatching = true);

      equation

        products.definition = productsSubstanceData;
        if nS>0 then
          products.solution = fill(substrates[1].solution,nP);
        end if;

        annotation (
          Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
                100,100}})),
          Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",       info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
<p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
<td><p>dissociation constant</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>G = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>G<sub>i</sub>) = &Delta;<sub>r</sub>H - T&middot;&Delta;<sub>r</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K) </p></td>
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
<h4><span style=\"color:#008000\">Notations</span></h4>
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
      end PartialReactionWithSubstanceData;

      partial model PartialReaction "Chemical Reaction"
        import Chemical;
        import Chemical.Utilities.Types.InitializationMethods;

        replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby Chemical.Interfaces.StateOfMatter
        "Substance model to translate data into substance properties"
          annotation (choices(
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.Incompressible  "Incompressible"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGas        "Ideal Gas"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
            choice(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

        parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
          annotation(Dialog(tab="Advanced"));
        parameter InitializationMethods initN_flow =Chemical.Utilities.Types.InitializationMethods.none  "Initialization method for n_flow"
          annotation(Dialog(tab= "Initialization", group="Molar flow"));
        parameter Modelica.Units.SI.MolarFlowRate n_flow_0 = 0 "Initial value for n_flow"
          annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.state)));
        parameter Utilities.Units.MolarFlowAcceleration n_acceleration_0 = 0 "Initial value for der(n_flow)"
          annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.derivative)));

        parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
        parameter Utilities.Units.Inertance L = dropOfCommons.L "Inertance of the flow"
          annotation(Dialog(tab="Advanced"));

        parameter Integer nS=0 "Number of substrate types"
          annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));

        parameter Modelica.Units.SI.StoichiometricNumber s[nS]=ones(nS)
          "Stoichiometric reaction coefficient for substrates"
          annotation (HideResult=true);

        parameter Integer nP=0 "Number of product types"
          annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));

        parameter Modelica.Units.SI.StoichiometricNumber p[nP]=ones(nP)
          "Stoichiometric reaction coefficients for products"
          annotation (HideResult=true);

        Modelica.Units.SI.MolarFlowRate rr(stateSelect=n_flowStateSelect) "Reaction molar flow rate";

        Chemical.Undirected.Interfaces.Rear substrates[nS](redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-100,0}), iconTransformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-100,0})));

        Chemical.Undirected.Interfaces.Fore products[nP](redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={100,0}), iconTransformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={100,0})));

        Modelica.Units.SI.MolarEnthalpy h_fore_mix, h_rear_mix;

        Real duRT_fore, duRT_rear, du_fore, du_rear, dr, Sx_fore,Px_rear,Kx;


        Modelica.Units.SI.ChemicalPotential uPure_substrates[nS];
        Modelica.Units.SI.ChemicalPotential uPure_products[nP];
      protected
        outer DropOfCommons dropOfCommons;

      public
          //debug only:
        Real Px_fore,Sx_rear;

      initial equation
        if initN_flow == InitializationMethods.state then
          rr = n_flow_0;
        elseif initN_flow == InitializationMethods.derivative then
          der(rr) = n_acceleration_0;
        elseif initN_flow == InitializationMethods.steadyState then
          der(rr) = 0;
        end if;

      equation

        du_fore = (s * substrates.state_forwards.u) - (p * products.state_forwards.u);
        du_rear = (s * substrates.state_rearwards.u) - (p * products.state_rearwards.u);

        duRT_fore = ((s * (substrates.state_forwards.u ./ (Modelica.Constants.R*substrates.solution.T))) - (p * (products.state_forwards.u ./ (Modelica.Constants.R*products.solution.T))));
        duRT_rear = ((s * (substrates.state_rearwards.u ./ (Modelica.Constants.R*substrates.solution.T))) - (p * (products.state_rearwards.u ./ (Modelica.Constants.R*products.solution.T))));

        for i in 1:nS loop
         uPure_substrates[i] = stateOfMatter.electroChemicalPotentialPure(
          substrates[i].definition,
          substrates[i].solution.T,
          substrates[i].solution.p,
          substrates[i].solution.v,
          substrates[i].solution.I);
        end for;

        for i in 1:nP loop
         uPure_products[i] = stateOfMatter.electroChemicalPotentialPure(
         products[i].definition,
         products[i].solution.T,
         products[i].solution.p,
         products[i].solution.v,
         products[i].solution.I);
        end for;

        Sx_fore = exp(s * ((substrates.state_forwards.u - uPure_substrates)./(Modelica.Constants.R*substrates.solution.T)));
        Px_rear = exp((p * ((products.state_rearwards.u - uPure_products)./(Modelica.Constants.R*products.solution.T))));

        //debug
        Sx_rear = exp(s * ((substrates.state_rearwards.u - uPure_substrates)./(Modelica.Constants.R*substrates.solution.T)));
        Px_fore = exp((p * ((products.state_forwards.u - uPure_products)./(Modelica.Constants.R*products.solution.T))));
        Kx = exp(- ((s * ((uPure_substrates)./(Modelica.Constants.R*substrates.solution.T))) - (p * ((uPure_products)./(Modelica.Constants.R*products.solution.T)))));

        //reaction molar rates
        rr*s = substrates.n_flow;
        rr*p = -products.n_flow;

        products.state_forwards.h = h_fore_mix*ones(nP);
        substrates.state_rearwards.h = h_rear_mix*ones(nS);


        if nS>0 then
          h_rear_mix*(substrates.n_flow*ones(nS)) + products.n_flow*products.state_rearwards.h = 0;
        else
          h_rear_mix = 0;
        end if;

        if nP>0 then
          h_fore_mix*(products.n_flow*ones(nP)) + substrates.n_flow*substrates.state_forwards.h = 0;
        else
          h_fore_mix = 0;
        end if;

        dr = (s * substrates.r) - (p * products.r);

        if nP>0 then

          (p * products.r) = (s * substrates.r)  -  der(rr)*L;

          for i in 2:nP loop
            //first product is based on inertial potential,
            //other products are provided as source with fixed flow and adaptation of their potential
            der(products[i].state_forwards.u).*TC = products[i].r;
          end for;
          for i in 2:nS loop
            der(substrates[i].state_rearwards.u).*TC = substrates[i].r;
          end for;
        end if;

        annotation (
          Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
                100,100}})),
          Documentation(revisions="<html>
<p><i>2013-2025 by </i>Marek Mateják </p>
</html>",       info="<html>
<h4><span style=\"color: #008000\">Notations</span></h4>
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
      end PartialReaction;
    end Internal;
    annotation (Documentation(info="<html>
<u>This package contains models implementing undirected versions of the processes. Here, the thermodynamic state of one or more fluid streams is changed by exchanging heat or work with the streams, or by delaying the state.</u>
</html>",   revisions="<html>
<u><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</u>
</html>"),   Icon(graphics={
           Ellipse(
            extent={{-60,54},{60,-66}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{-94,0},{94,0}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-64,60},{56,-60}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end Processes;

  package HeatExchangers
    extends Modelica.Icons.Package;

    model DiscretizedCrossFlowHEX "Discretized heat exchanger for single- or two-phase working fluid without potential drop"
      extends Internal.PartialDiscretizedHEX(nCellsParallel=nCells,crossFlow=true);

      Chemical.Undirected.Processes.FlowResistance flowResistanceA[nCells](
        redeclare package Medium = MediumA,
        each r(each displayUnit="mm") = 0.025,
        each l=1,
        redeclare function pLoss =
            Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
            (                                                                                                       each k=50))
        annotation (Placement(transformation(extent={{-20,-90},{-40,-70}})));
      Chemical.Undirected.Topology.JunctionMN junctionMN(
        redeclare package Medium = MediumA,
        N=1,
        M=nCells) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=0,
            origin={50,-80})));
      Chemical.Undirected.Topology.JunctionMN junctionMN1(redeclare package Medium = MediumA, N=nCells)
        annotation (Placement(transformation(extent={{-60,-90},{-80,-70}})));

    initial equation

      if initializeMassFlow then
        rearB.n_flow = n_flow_0_B;
        flowResistanceA.n_flow = n_flow_0_A/nCells*ones(nCells);
      else
        for i in 1:nCells - 1 loop
          flowResistanceA[i + 1].n_flow = flowResistanceA[1].n_flow;
        end for;
      end if;

    equation
      //Connecting equations (to interconnect pipes)

      //Fluid side B
      connect(rearB, thermalElementB[1].rear) annotation (Line(points={{-100,80},{-62,80},{-62,80},{-10,80}},
                                                                                            color={158,66,200}));
      for i in 1:nCells - 1 loop
        connect(thermalElementB[i].fore, thermalElementB[i + 1].rear);
      end for;
      connect(thermalElementB[nCells].fore, foreB) annotation (Line(points={{10,80},{62,80},{62,80},{100,80}},
                                                                                                   color={158,66,200}));

      connect(thermalElementA.heatPort, thermalConductor.port_a) annotation (Line(points={{4.44089e-16,-70.2},{4.44089e-16,-40},{0,-40},{0,-10}}, color={191,0,0}));
      connect(thermalElementB.heatPort, thermalConductor.port_b) annotation (Line(points={{4.44089e-16,70.2},{4.44089e-16,40},{0,40},{0,10}}, color={191,0,0}));

      //Fluid side A
      connect(junctionMN1.fores[1], foreA) annotation (Line(
          points={{-80,-80},{-100,-80}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.rears[1], rearA) annotation (Line(
          points={{60,-80},{100,-80}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.fores, thermalElementA.rear) annotation (Line(
          points={{40,-80},{10,-80}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistanceA.rear, thermalElementA.fore) annotation (Line(
          points={{-20,-80},{-10,-80}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistanceA.fore, junctionMN1.rears) annotation (Line(
          points={{-40,-80},{-60,-80}},
          color={158,66,200},
          thickness=0.5));

      annotation (Icon(graphics={
            Text(
              extent={{-72,76},{-60,64}},
              textColor={158,66,200},
              textString="N"),
            Text(
              extent={{-42,76},{-30,64}},
              textColor={158,66,200},
              textString="..."),
            Text(
              extent={{-10,76},{2,64}},
              textColor={158,66,200},
              textString="..."),
            Text(
              extent={{20,76},{32,64}},
              textColor={158,66,200},
              textString="2"),
            Text(
              extent={{50,76},{62,64}},
              textColor={158,66,200},
              textString="1"),
            Text(
              extent={{10,134},{50,94}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{80,-94},{120,-134}},
              textColor={175,175,175},
              textString="B")}), Documentation(info="<html>
<u>The undirected cross-flow discretized heat exchanger uses a number of conduction elements (which is set by the parameter nCells) as discrete control volumes to exchange heat between two fluid streams. </u>
<u>Side A splits the fluid stream into nCells substreams that are parallel. The flow-resistance is chosen to be very small and only ensures numerical stability of the parallel streams. By default, it is a linear-quadratic flow resistance, so the mass flow through each of the parallel streams is the same. If exchanged for flow-resistance that depends on media properties (e.g. a laminar-turbulent) the mass flow on the paths will be different. For side B the elements are serial and numbered 1 to nCells in the flow direction. The elements&apos; heatports are connected via a thermal conductor that models the wall. The connections are ordered to result in a cross-flow configuration. </u>
<u>The conduction elements are computing a heat transfer coefficient between their heatport and the fluid contained. They are replaceable with a choice between a single-phase and a two-phase version, both can be further parametrized. Although the single-phase version works for two-phase media (not the other way around), using the two-phase one for two-phase media enables to set different heat transfer coefficients depending on the phase (liquid/gaseous/2-phase) state of the medium. </u>
<u>Note that since the model uses conductionElements as discrete control volumes that in turn assume quasi-stationary mass and therefore are part of a fluid stream rather than break it into two (like a full volume would), the same holds for both sides of the heat exchanger - they are part of a fluid stream and don&apos;t break it. The quasi-stationary mass assumption also implies that for (fast) changing masses/densities in any of the conduction elements the heat exchanger will (slightly) violate the conservation of energy. Furthermore, the conduction elements change their behavior for reversed mass flow, therefore this model asserts for negative mass flow with the level dropOfCommons.assertionLevel. </u>
<u>The parameters A (heat transferring area), k_wall (heat transfer coefficient of the wall between the streams) and the heat transfer coefficients in the conduction elements scale the transferred heat (the middle only if the wall and the latter only of the heat transfer into a fluid is the choke of the heatflow). </u>
<u>The parameter V determines the amount of fluid in the heat exchanger and therefore the dynamic for non-steady states. </u>
<u>The initialization tab allows for a mass flow initialization for both paths, as well as to determine from which direction the enthalpy in the control volumes should be initialized (fore/rear), or if it should start with a given enthalpy. The other option is to initialize the enthalpy with a given value. </u>
<u>The Advanced tab allows to modify the mass flow that triggers the reverse-mass-flow-assertion and has an option to enforce global conservation of energy. The latter is done by feeding back any energy the conduction elements accumulated over time, basically making it impossible to store energy in their fluid long-term. While this enforces long-term conservation of energy it changes the medium-/short-term dynamics of the system and is therefore disabled by default. </u>
</html>"));
    end DiscretizedCrossFlowHEX;

    model DiscretizedCounterFlowHEX "Discretized heat exchanger for single- or two-phase working fluids without potential drop"
      extends Internal.PartialDiscretizedHEX;

    initial equation

      if initializeMassFlow then
        rearA.n_flow = n_flow_0_A;
        rearB.n_flow = n_flow_0_B;
      end if;

    equation

      //Connecting equations (to interconnect pipes)
      //Fluid side B
      connect(rearB, thermalElementB[1].rear) annotation (Line(points={{-100,80},{-10,80}}, color={158,66,200}));
      for i in 1:nCells - 1 loop
        connect(thermalElementB[i].fore, thermalElementB[i + 1].rear);
      end for;
      connect(thermalElementB[nCells].fore, foreB) annotation (Line(points={{10,80},{100,80}}, color={158,66,200}));

      //Fluid side A
      connect(rearA, thermalElementA[1].rear) annotation (Line(points={{100,-80},{10,-80}}, color={158,66,200}));
      for i in 1:nCells - 1 loop
        connect(thermalElementA[i].fore, thermalElementA[i + 1].rear);
      end for;
      connect(thermalElementA[nCells].fore, foreA) annotation (Line(points={{-10,-80},{-100,-80}}, color={158,66,200}));

      connect(thermalElementB.heatPort, thermalConductor.port_b) annotation (Line(points={{0,70.2},{0,10}}, color={191,0,0}));

      for i in 1:nCells loop
        connect(thermalElementA[i].heatPort, thermalConductor[nCells + 1 - i].port_a) annotation (Line(points={{-6.66134e-16,-70.2},{-6.66134e-16,-10},{0,-10}}, color={191,0,0}));
      end for;

      annotation (Icon(graphics={
            Text(
              extent={{-70,76},{-58,64}},
              textColor={158,66,200},
              textString="1"),
            Text(
              extent={{-40,76},{-28,64}},
              textColor={158,66,200},
              textString="2"),
            Text(
              extent={{-8,76},{4,64}},
              textColor={158,66,200},
              textString="..."),
            Text(
              extent={{22,76},{34,64}},
              textColor={158,66,200},
              textString="..."),
            Text(
              extent={{50,76},{62,64}},
              textColor={158,66,200},
              textString="N"),
            Text(
              extent={{-120,132},{-80,92}},
              textColor={175,175,175},
              textString="B"),
            Text(
              extent={{80,-94},{120,-134}},
              textColor={175,175,175},
              textString="A")}), Documentation(info="<html>
<u>The undirected counter-flow discretized heat exchanger uses a number of conduction elements (which is set by the parameter nCells) as discrete control volumes to exchange heat between two fluid streams. </u>
<u>For each side the elements are numbered 1 to nCells from rear to fore and the elements&apos; heatports are connected via a thermal conductor that models the wall. The connections are ordered to result in a counter-flow configuration. </u>
<u>The conduction elements are computing a heat transfer coefficient between their heatport and the fluid contained. They are replaceable with a choice between a single-phase and a two-phase version, both can be further parametrized. Although the single-phase version works for two-phase media (not the other way around), using the two-phase one for two-phase media enables to set different heat transfer coefficients depending on the phase (liquid/gaseous/2-phase) state of the medium. </u>
<u>Note that since the model uses conductionElements as discrete control volumes that in turn assume quasi-stationary mass and therefore are part of a fluid stream rather than break it into two (like a full volume would), the same holds for both sides of the heat exchanger; they are part of a fluid stream and don&apos;t break it. The quasi-stationary mass assumption also implies that for (fast) changing masses/densities in any of the conduction elements the heat exchanger will (slightly) violate the conservation of energy.</u>
<u>The parameters A (heat transferring area), k_wall (heat transfer coefficient of the wall between the streams) and the heat transfer coefficients in the conduction elements scale the transferred heat (the middle only if the wall and the latter only of the heat transfer into a fluid is the choke of the heatflow). </u>
<u>The parameter V determines the amount of fluid in the heat exchanger and therefore the dynamic for non-steady states. </u>
<u>The initialization tab allows for a mass flow initialization for both paths, as well as to determine from which direction the enthalpy in the control volumes should be initialized (fore/rear), or if it should start with a given enthalpy. The other option is to initialize the enthalpy with a given value.</u>
<u>The Advanced tab allows to influence the mass flow regularization for near zero mass flow and has an option to enforce global conservation of energy. The latter is done by feeding back any energy the conduction elements accumulated over time, basically making it impossible to store energy in their fluid long-term. While this enforces long-term conservation of energy it changes the medium-/short-term dynamics of the system and is therefore disabled by default. </u>
</html>"));
    end DiscretizedCounterFlowHEX;

    package Tests
      extends Modelica.Icons.ExamplesPackage;

      model TestDiscretizedHEX

        extends Modelica.Icons.Example;

        replaceable package MediumAir = Chemical.Media.myMedia.Air.MoistAir constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium model" annotation (
            choicesAllMatching=true,
            Dialog(group = "Medium definitions"));

        replaceable package MediumRefrigerant =
            Chemical.Media.myMedia.R134a.R134a_ph                                     constrainedby
          Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium model" annotation (
            choicesAllMatching=true,
            Dialog(group = "Medium definitions"));

        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package Medium = MediumAir,
          potentialFromInput=true,
          T0_par=311.15) annotation (Placement(transformation(extent={{-96,14},{-76,34}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(
          redeclare package Medium = MediumAir,
          potentialFromInput=true,
          T0_par=311.15,
          u0_par=100000) annotation (Placement(transformation(extent={{84,14},{104,34}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore1(
          redeclare package Medium = MediumRefrigerant,
          setEnthalpy=true,
          potentialFromInput=true,
          h0_par=300e3,
          T0_par=278.15,
          u0_par=400000) annotation (Placement(transformation(extent={{-90,-50},{-110,-30}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear1(
          redeclare package Medium = MediumRefrigerant,
          setEnthalpy=true,
          temperatureFromInput=false,
          potentialFromInput=true,
          h0_par=450e3,
          T0_par=268.15) annotation (Placement(transformation(extent={{112,-18},{92,2}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm(
          redeclare package Medium = MediumAir,
          temperatureUnit="degC",
          potentialUnit="bar") annotation (Placement(transformation(extent={{-42,22},{-22,42}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm1(
          redeclare package Medium = MediumAir,
          temperatureUnit="degC",
          potentialUnit="bar",
          outputMassFlowRate=true) annotation (Placement(transformation(extent={{18,22},{38,42}})));
        Chemical.Undirected.Processes.FlowResistance flowResistanceA(
          redeclare package Medium = MediumAir,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=1,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-70,14},{-50,34}})));
        inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{74,74},{94,94}})));
        Modelica.Blocks.Sources.Ramp ramp1(
          height=-0.6,
          duration=10,
          offset=0.3,
          startTime=15) annotation (Placement(transformation(extent={{-40,-58},{-20,-38}})));
        Modelica.Blocks.Sources.Constant const1(k=1e5) annotation (Placement(transformation(extent={{-120,14},{-100,34}})));
        DiscretizedCounterFlowHEX discretizedHEX(
          redeclare package MediumA = MediumAir,
          redeclare package MediumB = MediumRefrigerant,
          redeclare model ConductionElementA = Internal.ConductionElementHEX,
          redeclare model ConductionElementB =
              Internal.ConductionElementHEX_twoPhase,
          nCells=10,
          V_Hex(displayUnit="m3"),
          initializeMassFlow=false,
          k_wall=300) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-4,16})));
        Chemical.Undirected.Processes.FlowResistance flowResistanceB(
          redeclare package Medium = MediumRefrigerant,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=0.3,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-52,-2},{-72,18}})));
        Modelica.Blocks.Continuous.PI PI1(
          k=-10000,
          T=0.1,
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          y_start=1e5)
          annotation (Placement(transformation(extent={{94,46},{114,66}})));
        Modelica.Blocks.Nonlinear.Limiter limiter1(uMax=5e5, uMin=100)
          annotation (Placement(transformation(extent={{126,18},{114,30}})));
        Modelica.Blocks.Math.Feedback feedback1
          annotation (Placement(transformation(extent={{38,46},{58,66}})));
        Modelica.Blocks.Sources.RealExpression airFlow_setPoint1(y=1)
          annotation (Placement(transformation(extent={{6,72},{26,92}})));
        Modelica.Blocks.Sources.Ramp ramp2(
          height=-1,
          duration=1,
          offset=1,
          startTime=30) annotation (Placement(transformation(extent={{-22,54},{-2,74}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm2(
          redeclare package Medium = MediumRefrigerant,
          temperatureUnit="degC",
          potentialUnit="bar",
          outputMassFlowRate=true) annotation (Placement(transformation(extent={{38,10},{18,-10}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm3(
          redeclare package Medium = MediumRefrigerant,
          temperatureUnit="degC",
          potentialUnit="bar") annotation (Placement(transformation(extent={{-24,10},{-44,-10}})));
        Modelica.Blocks.Continuous.PI PI(
          k=10000,
          T=0.001,
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          x_start=300,
          y_start=30e5)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={50,-48})));
        Modelica.Blocks.Math.Feedback feedback
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=0,
              origin={8,-48})));
        Modelica.Blocks.Nonlinear.Limiter limiter(uMax=35e5, uMin=1e5)
          annotation (Placement(transformation(extent={{-6,-6},{6,6}},
              rotation=0,
              origin={84,-48})));
        Modelica.Blocks.Sources.Ramp ramp3(
          height=-26e5,
          duration=10,
          offset=30e5,
          startTime=15) annotation (Placement(transformation(extent={{-134,-50},{-114,-30}})));
      equation
        connect(boundary_rear.u0_var, const1.y) annotation (Line(points={{-88,30},{-94,30},{-94,24},{-99,24}},
                                                                                             color={0,0,127}));
        connect(boundary_fore1.rear, flowResistanceB.fore) annotation (Line(
            points={{-90,-40},{-84,-40},{-84,8},{-72,8}},
            color={158,66,200},
            thickness=0.5));
        connect(PI1.y,limiter1. u) annotation (Line(points={{115,56},{136,56},{136,24},{127.2,24}},
                                   color={0,0,127}));
        connect(boundary_fore.u0_var, limiter1.y) annotation (Line(points={{96,30},{104,30},{104,24},{113.4,24}},
                                                                                                color={0,0,127}));
        connect(multiSensor_Tpm1.n_flow_out, feedback1.u2) annotation (Line(points={{38,28},{48,28},{48,48}}, color={0,0,127}));
        connect(feedback1.y, PI1.u) annotation (Line(points={{57,56},{92,56}}, color={0,0,127}));
        connect(boundary_rear1.fore, multiSensor_Tpm2.rear)
          annotation (Line(
            points={{92,-8},{66,-8},{66,8},{38,8}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm3.fore, flowResistanceB.rear) annotation (Line(
            points={{-44,8},{-52,8}},
            color={158,66,200},
            thickness=0.5));
        connect(feedback.y,PI. u)
          annotation (Line(points={{17,-48},{38,-48}}, color={0,0,127}));
        connect(PI.y,limiter. u)
          annotation (Line(points={{61,-48},{76.8,-48}}, color={0,0,127}));
        connect(limiter.y, boundary_rear1.u0_var) annotation (Line(points={{90.6,-48},{112,-48},{112,-2},{104,-2}}, color={0,0,127}));
        connect(multiSensor_Tpm2.n_flow_out, feedback.u2) annotation (Line(points={{18,4},{8,4},{8,-40}}, color={0,0,127}));
        connect(multiSensor_Tpm1.fore, boundary_fore.rear) annotation (Line(
            points={{38,24},{84,24}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear.fore, flowResistanceA.rear) annotation (Line(
            points={{-76,24},{-70,24}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm.rear, flowResistanceA.fore) annotation (Line(
            points={{-42,24},{-50,24}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore1.u0_var, ramp3.y) annotation (Line(points={{-102,-34},{-108,-34},{-108,-40},{-113,-40}}, color={0,0,127}));
        connect(ramp2.y, feedback1.u1) annotation (Line(points={{-1,64},{24,64},{24,56},{40,56}}, color={0,0,127}));
        connect(feedback.u1, ramp1.y) annotation (Line(points={{0,-48},{-19,-48}}, color={0,0,127}));
        connect(discretizedHEX.foreB, multiSensor_Tpm3.rear) annotation (Line(
            points={{-14.2,8},{-24,8}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm.fore, discretizedHEX.rearA) annotation (Line(
            points={{-22,24},{-14,24}},
            color={158,66,200},
            thickness=0.5));
        connect(discretizedHEX.foreA, multiSensor_Tpm1.rear) annotation (Line(
            points={{6,24},{18,24}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm2.fore, discretizedHEX.rearB) annotation (Line(
            points={{18,8},{6,8}},
            color={158,66,200},
            thickness=0.5));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=40, Tolerance=1e-6, Interval=0.04, __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
        <u>Owner: <a href=\"mailto:niels.weber@dlr.de\">Niels Weber</a></u>
</html>"));
      end TestDiscretizedHEX;

      model ConductionElementTwoPhase
        extends Modelica.Icons.Example;

        package MediumRefrigerant = Chemical.Media.myMedia.R134a.R134a_ph;

        Internal.ConductionElementHEX_twoPhase conductionElementHEX_twoPhase(
          redeclare package Medium = MediumRefrigerant,
          V(displayUnit="l") = 0.001,
          init=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.h,
          A=10,
          U_liq_nom=700,
          U_vau_nom=500,
          U_tu_nom=1000) annotation (Placement(transformation(extent={{-12,-10},{8,10}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package Medium = MediumRefrigerant,
          setEnthalpy=true,
          enthalpyFromInput=false,
          h0_par=200e3,
          u0_par=500000) annotation (Placement(transformation(extent={{-82,-10},{-62,10}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(
          redeclare package Medium = MediumRefrigerant,
          setEnthalpy=true,
          potentialFromInput=true,
          h0_par=450e3) annotation (Placement(transformation(extent={{60,-10},{80,10}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance(
          redeclare package Medium = MediumRefrigerant,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-44,-10},{-24,10}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm(
          redeclare package Medium = MediumRefrigerant,
          temperatureUnit="degC",
          outputMassFlowRate=true) annotation (Placement(transformation(extent={{24,-2},{44,18}})));
        Modelica.Blocks.Continuous.PI PI(
          k=-10000,
          T=0.1,
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          x_start=300,
          y_start=5e5)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={90,34})));
        Modelica.Blocks.Math.Feedback feedback
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={50,46})));
        Modelica.Blocks.Nonlinear.Limiter limiter(uMax=35e5, uMin=1e5)
          annotation (Placement(transformation(extent={{6,-6},{-6,6}},
              rotation=0,
              origin={88,0})));
        Modelica.Blocks.Sources.Ramp ramp1(
          height=-0.6,
          duration=5,
          offset=0.3,
          startTime=20) annotation (Placement(transformation(extent={{-4,36},{16,56}})));
      equation
        connect(boundary_rear.fore, flowResistance.rear) annotation (Line(
            points={{-62,0},{-44,0}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance.fore, conductionElementHEX_twoPhase.rear)
          annotation (Line(
            points={{-24,0},{-12,0}},
            color={158,66,200},
            thickness=0.5));
        connect(conductionElementHEX_twoPhase.fore, multiSensor_Tpm.rear)
          annotation (Line(
            points={{8,0},{24,0}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore.rear, multiSensor_Tpm.fore) annotation (Line(
            points={{60,0},{44,0}},
            color={158,66,200},
            thickness=0.5));
        connect(feedback.y,PI. u)
          annotation (Line(points={{59,46},{68,46},{68,34},{78,34}},
                                                        color={0,0,127}));
        connect(PI.y,limiter. u)
          annotation (Line(points={{101,34},{110,34},{110,0},{95.2,0}},
                                                         color={0,0,127}));
        connect(boundary_fore.u0_var, limiter.y) annotation (Line(points={{72,6},{76,6},{76,0},{81.4,0}},
                                                                                            color={0,0,127}));
        connect(multiSensor_Tpm.n_flow_out, feedback.u2) annotation (Line(points={{44,4},{50,4},{50,38}}, color={0,0,127}));
        connect(ramp1.y, feedback.u1) annotation (Line(points={{17,46},{42,46}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=100, Tolerance=1e-6, Interval=0.1, __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
        <u>Owner: <a href=\"mailto:niels.weber@dlr.de\">Niels Weber</a></u>
</html>"));
      end ConductionElementTwoPhase;

      model TestDiscretizedHEXvsDir

        extends Modelica.Icons.Example;

        replaceable package MediumAir = Chemical.Media.myMedia.Air.DryAirNasa constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium model" annotation (
            choicesAllMatching=true,
            Dialog(group = "Medium definitions"));

        replaceable package MediumRefrigerant =
            Chemical.Media.myMedia.R134a.R134a_ph                                     constrainedby
          Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium model" annotation (
            choicesAllMatching=true,
            Dialog(group = "Medium definitions"));

        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package Medium = MediumAir,
          potentialFromInput=true,
          T0_par=311.15) annotation (Placement(transformation(extent={{-94,198},{-74,218}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(
          redeclare package Medium = MediumAir,
          potentialFromInput=true,
          T0_par=303.15,
          u0_par=100000) annotation (Placement(transformation(extent={{86,198},{106,218}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore1(
          redeclare package Medium = MediumRefrigerant,
          setEnthalpy=true,
          potentialFromInput=true,
          h0_par=300e3,
          T0_par=278.15,
          u0_par=400000) annotation (Placement(transformation(extent={{-88,134},{-108,154}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear1(
          redeclare package Medium = MediumRefrigerant,
          setEnthalpy=true,
          temperatureFromInput=false,
          potentialFromInput=true,
          h0_par=450e3,
          T0_par=268.15) annotation (Placement(transformation(extent={{114,166},{94,186}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm(
          redeclare package Medium = MediumAir,
          temperatureUnit="degC",
          potentialUnit="bar") annotation (Placement(transformation(extent={{-40,206},{-20,226}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm1(
          redeclare package Medium = MediumAir,
          temperatureUnit="degC",
          potentialUnit="bar",
          outputMassFlowRate=true) annotation (Placement(transformation(extent={{20,206},{40,226}})));
        Chemical.Undirected.Processes.FlowResistance flowResistanceA(
          redeclare package Medium = MediumAir,
          n_flow_0=0,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-68,198},{-48,218}})));
        Modelica.Blocks.Sources.Ramp rampChemicalPotential(
          height=1e5,
          duration=1,
          offset=1e5,
          startTime=15) annotation (Placement(transformation(extent={{-134,186},{-114,206}})));
        inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{76,258},{96,278}})));
        Modelica.Blocks.Sources.Ramp ramp1(
          height=-0.5,
          duration=1,
          offset=0.3,
          startTime=15) annotation (Placement(transformation(extent={{-34,126},{-14,146}})));
        DiscretizedCounterFlowHEX discretizedHEXUndir(
          redeclare package MediumA = MediumAir,
          redeclare package MediumB = MediumRefrigerant,
          redeclare model ConductionElementA = Internal.ConductionElementHEX,
          redeclare model ConductionElementB =
              Internal.ConductionElementHEX_twoPhase,
          nCells=10,
          V_Hex(displayUnit="m3"),
          initializeMassFlow=true,
          k_wall=300,
          calculate_efficiency=true)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-2,200})));
        Chemical.Undirected.Processes.FlowResistance flowResistanceB(
          redeclare package Medium = MediumRefrigerant,
          n_flow_0=0,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-50,182},{-70,202}})));
        Modelica.Blocks.Continuous.PI PI1(
          k=-10000,
          T=0.1,
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          y_start=1e5)
          annotation (Placement(transformation(extent={{96,230},{116,250}})));
        Modelica.Blocks.Nonlinear.Limiter limiter1(uMax=5e5, uMin=100)
          annotation (Placement(transformation(extent={{128,202},{116,214}})));
        Modelica.Blocks.Math.Feedback feedback1
          annotation (Placement(transformation(extent={{40,230},{60,250}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm2(
          redeclare package Medium = MediumRefrigerant,
          temperatureUnit="degC",
          outputMassFlowRate=true) annotation (Placement(transformation(extent={{40,194},{20,174}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm3(redeclare
            package                                                                    Medium =
              MediumRefrigerant,                                                                                   temperatureUnit="degC")
          annotation (Placement(transformation(extent={{-22,194},{-42,174}})));
        Modelica.Blocks.Continuous.PI PI(
          k=10000,
          T=0.001,
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          x_start=300,
          y_start=30e5)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={52,136})));
        Modelica.Blocks.Math.Feedback feedback
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=0,
              origin={10,136})));
        Modelica.Blocks.Nonlinear.Limiter limiter(uMax=35e5, uMin=1e5)
          annotation (Placement(transformation(extent={{-6,-6},{6,6}},
              rotation=0,
              origin={86,136})));
        Modelica.Blocks.Sources.Ramp ramp3(
          height=-26e5,
          duration=1,
          offset=30e5,
          startTime=15) annotation (Placement(transformation(extent={{-144,134},{-124,154}})));
        Chemical.Boundaries.Source sourceA(
          redeclare package Medium = MediumAir,
          T0_par=303.15,
          u0_par=200000)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={112,-64})));
        Chemical.Boundaries.Sink sinkA(
          redeclare package Medium = MediumAir,
          potentialFromInput=true,
          u0_par=100000) annotation (Placement(transformation(extent={{-44,-74},{-64,-54}})));
        Chemical.Sensors.MultiSensor_Tpm multiSensor_Tpm4(
          redeclare package Medium = MediumAir,
          temperatureUnit="degC",
          potentialUnit="bar") annotation (
            Placement(transformation(
              extent={{11,10},{-11,-10}},
              rotation=0,
              origin={49,-74})));
        Chemical.Sensors.MultiSensor_Tpm multiSensor_Tpm5(
          redeclare package Medium = MediumAir,
          digits=3,
          outputMassFlowRate=true,
          temperatureUnit="degC") annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-22,-74})));
        Chemical.Boundaries.Source sourceB(
          redeclare package Medium = MediumRefrigerant,
          setEnthalpy=true,
          temperatureFromInput=false,
          potentialFromInput=true,
          T0_par=283.15,
          u0_par=200000,
          h0_par=300e3)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-58,-20})));
        Chemical.Boundaries.Sink sinkB(
          redeclare package Medium = MediumRefrigerant,
          potentialFromInput=false,
          u0_par(displayUnit="bar") = 400000)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=90,
              origin={64,2})));
        Chemical.Sensors.MultiSensor_Tpm multiSensor_Tpm6(
          redeclare package Medium = MediumRefrigerant,
          digits=3,
          temperatureUnit="degC")
          annotation (Placement(transformation(extent={{26,-48},{46,-28}})));
        Modelica.Blocks.Sources.RealExpression refFlow_setPoint1(y=0.2)
          annotation (Placement(transformation(extent={{10,10},{-10,-10}},
              rotation=0,
              origin={26,16})));
        Chemical.Sensors.MultiSensor_Tpm multiSensor_Tpm7(
          redeclare package Medium = MediumRefrigerant,
          digits=3,
          outputMassFlowRate=true,
          temperatureUnit="degC")
          annotation (Placement(transformation(extent={{-38,-48},{-18,-28}})));
        Modelica.Blocks.Sources.RealExpression airFlow_setPoint2(y=1)
          annotation (Placement(transformation(extent={{-160,-58},{-140,-38}})));
        Modelica.Blocks.Continuous.PI PI2(
          k=10000,
          T=0.1,
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          x_start=40,
          y_start=4e5)
          annotation (Placement(transformation(extent={{10,10},{-10,-10}},
              rotation=0,
              origin={-32,16})));
        Modelica.Blocks.Math.Feedback feedback2
          annotation (Placement(transformation(extent={{10,-10},{-10,10}},
              rotation=0,
              origin={-4,16})));
        Modelica.Blocks.Nonlinear.Limiter limiter2(uMax=10e5, uMin=1e5)
          annotation (Placement(transformation(extent={{6,6},{-6,-6}},
              rotation=90,
              origin={-58,-2})));
        Modelica.Blocks.Continuous.PI PI3(
          k=-10000,
          T=0.1,
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          y_start=1e5)
          annotation (Placement(transformation(extent={{-106,-74},{-86,-54}})));
        Modelica.Blocks.Math.Feedback feedback3
          annotation (Placement(transformation(extent={{-134,-74},{-114,-54}})));
        Modelica.Blocks.Nonlinear.Limiter limiter3(uMax=5e5, uMin=100)
          annotation (Placement(transformation(extent={{-76,-70},{-64,-58}})));
        Chemical.HeatExchangers.DiscretizedCounterFlowHEX evaporator(
          redeclare model ConductionElementB =
              Chemical.HeatExchangers.Internal.ConductionElementHEX_twoPhase,
          redeclare package MediumA = MediumAir,
          redeclare package MediumB = MediumRefrigerant,
          initializeMassFlow=true,
          nCells=10,
          k_wall=300,
          calculate_efficiency=true)
          annotation (Placement(transformation(extent={{10,10},{-10,-10}},
              rotation=180,
              origin={2,-56})));
        Chemical.Processes.FlowResistance flowResistanceA1(
          redeclare package Medium = MediumAir,
          n_flow_0=1,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (
            material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{90,-74},{70,-54}})));
        Chemical.Processes.FlowResistance flowResistanceB1(
          redeclare package Medium = MediumRefrigerant,
          n_flow_0=0.3,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (
            material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={64,-26})));
        Chemical.Sensors.SingleSensorSelect singleSensorSelect(
          redeclare package Medium = MediumRefrigerant,
          quantity=Chemical.Sensors.Internal.Types.Quantities.h_Jpkg)
          annotation (Placement(transformation(extent={{26,-14},{46,-34}})));
        Chemical.Sensors.SingleSensorSelect singleSensorSelect1(
          redeclare package Medium = MediumRefrigerant,
          quantity=Chemical.Sensors.Internal.Types.Quantities.h_Jpkg)
          annotation (Placement(transformation(extent={{-18,-14},{-38,-34}})));
        Chemical.Sensors.TwoPhaseSensorSelect sensorVaporQuality(
          redeclare package Medium = MediumRefrigerant,
          quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg)
          annotation (Placement(transformation(extent={{-18,-66},{-38,-46}})));
        Chemical.Sensors.TwoPhaseSensorSelect sensorVaporQuality1(
          redeclare package Medium = MediumRefrigerant,
          quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg)
          annotation (Placement(transformation(extent={{26,-66},{46,-46}})));
        Chemical.Boundaries.Source sourceA1(
          redeclare package Medium = MediumAir,
          T0_par=311.15,
          u0_par=100000)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-112,-166})));
        Chemical.Boundaries.Sink sinkA1(
          redeclare package Medium = MediumAir,
          potentialFromInput=true,
          u0_par=100000) annotation (Placement(transformation(extent={{102,-176},{122,-156}})));
        Chemical.Sensors.MultiSensor_Tpm multiSensor_Tpm8(
          redeclare package Medium = MediumAir,
          temperatureUnit="degC",
          potentialUnit="bar") annotation (
            Placement(transformation(
              extent={{-11,-10},{11,10}},
              rotation=0,
              origin={-37,-156})));
        Chemical.Sensors.MultiSensor_Tpm multiSensor_Tpm9(
          redeclare package Medium = MediumAir,
          digits=3,
          outputMassFlowRate=true,
          temperatureUnit="degC") annotation (Placement(transformation(
              extent={{10,10},{-10,-10}},
              rotation=180,
              origin={52,-156})));
        Chemical.Boundaries.Source sourceB1(
          redeclare package Medium = MediumRefrigerant,
          setEnthalpy=true,
          temperatureFromInput=false,
          potentialFromInput=true,
          T0_par=283.15,
          u0_par=200000,
          h0_par=450e3)
          annotation (Placement(transformation(extent={{10,-10},{-10,10}},
              rotation=270,
              origin={80,-206})));
        Chemical.Boundaries.Sink sinkB1(
          redeclare package Medium = MediumRefrigerant,
          potentialFromInput=false,
          u0_par(displayUnit="bar") = 3000000)
          annotation (Placement(transformation(extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-76,-228})));
        Chemical.Sensors.MultiSensor_Tpm multiSensor_Tpm10(
          redeclare package Medium = MediumRefrigerant,
          digits=3,
          temperatureUnit="degC")
          annotation (Placement(transformation(extent={{-24,-182},{-44,-202}})));
        Modelica.Blocks.Sources.RealExpression refFlow_setPoint2(y=0.3)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-14,-244})));
        Chemical.Sensors.MultiSensor_Tpm multiSensor_Tpm11(
          redeclare package Medium = MediumRefrigerant,
          digits=3,
          outputMassFlowRate=true,
          temperatureUnit="degC")
          annotation (Placement(transformation(extent={{40,-182},{20,-202}})));
        Modelica.Blocks.Sources.RealExpression airFlow_setPoint3(y=1)
          annotation (Placement(transformation(extent={{8,-140},{28,-120}})));
        Modelica.Blocks.Continuous.PI PI4(
          k=10000,
          T=0.001,
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          x_start=300,
          y_start=30e5)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={44,-244})));
        Modelica.Blocks.Math.Feedback feedback4
          annotation (Placement(transformation(extent={{-10,10},{10,-10}},
              rotation=0,
              origin={20,-244})));
        Modelica.Blocks.Nonlinear.Limiter limiter4(uMax=35e5, uMin=20e5)
          annotation (Placement(transformation(extent={{-6,-6},{6,6}},
              rotation=0,
              origin={68,-244})));
        Modelica.Blocks.Continuous.PI PI5(
          k=-10000,
          T=0.1,
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          y_start=1e5)
          annotation (Placement(transformation(extent={{108,-140},{128,-120}})));
        Modelica.Blocks.Math.Feedback feedback5
          annotation (Placement(transformation(extent={{68,-140},{88,-120}})));
        Modelica.Blocks.Nonlinear.Limiter limiter5(uMax=5e5, uMin=100)
          annotation (Placement(transformation(extent={{140,-172},{128,-160}})));
        Chemical.HeatExchangers.DiscretizedCounterFlowHEX condenser(
          redeclare package MediumA = MediumAir,
          redeclare package MediumB = MediumRefrigerant,
          redeclare model ConductionElementB =
              Chemical.HeatExchangers.Internal.ConductionElementHEX_twoPhase,
          initializeMassFlow=true,
          nCells=10,
          k_wall=300,
          calculate_efficiency=true)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-2,-174})));
        Chemical.Processes.FlowResistance flowResistanceA2(
          redeclare package Medium = MediumAir,
          n_flow_0=0.5,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (
            material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-84,-176},{-64,-156}})));
        Chemical.Processes.FlowResistance flowResistanceB2(
          redeclare package Medium = MediumRefrigerant,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.none,
          n_flow_0=0.3,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (
            material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-76,-202})));
        Chemical.Sensors.SingleSensorSelect singleSensorSelect2(
          redeclare package Medium = MediumRefrigerant,
          quantity=Chemical.Sensors.Internal.Types.Quantities.h_Jpkg)
          annotation (Placement(transformation(extent={{-24,-218},{-44,-198}})));
        Chemical.Sensors.SingleSensorSelect singleSensorSelect3(
          redeclare package Medium = MediumRefrigerant,
          quantity=Chemical.Sensors.Internal.Types.Quantities.h_Jpkg)
          annotation (Placement(transformation(extent={{22,-216},{42,-196}})));
        Chemical.Sensors.TwoPhaseSensorSelect sensorVaporQuality2(
          redeclare package Medium = MediumRefrigerant,
          quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg)
          annotation (Placement(transformation(extent={{-24,-230},{-44,-210}})));
        Chemical.Sensors.TwoPhaseSensorSelect sensorVaporQuality3(
          redeclare package Medium = MediumRefrigerant,
          quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg)
          annotation (Placement(transformation(extent={{22,-230},{42,-210}})));
        Modelica.Blocks.Sources.Ramp rampMassflow(
          height=-2,
          duration=1,
          offset=1,
          startTime=15) annotation (Placement(transformation(extent={{-6,234},{14,254}})));
      equation
        connect(boundary_fore1.rear, flowResistanceB.fore) annotation (Line(
            points={{-88,144},{-82,144},{-82,192},{-70,192}},
            color={158,66,200},
            thickness=0.5));
        connect(PI1.y,limiter1. u) annotation (Line(points={{117,240},{138,240},{138,208},{129.2,208}},
                                   color={0,0,127}));
        connect(boundary_fore.u0_var, limiter1.y) annotation (Line(points={{98,214},{106,214},{106,208},{115.4,208}},
                                                                                                color={0,0,127}));
        connect(multiSensor_Tpm1.n_flow_out, feedback1.u2) annotation (Line(points={{40,212},{50,212},{50,232}},
                                                                                                              color={0,0,127}));
        connect(feedback1.y, PI1.u) annotation (Line(points={{59,240},{94,240}},
                                                                               color={0,0,127}));
        connect(boundary_rear1.fore, multiSensor_Tpm2.rear)
          annotation (Line(
            points={{94,176},{68,176},{68,192},{40,192}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm3.fore, flowResistanceB.rear) annotation (Line(
            points={{-42,192},{-50,192}},
            color={158,66,200},
            thickness=0.5));
        connect(feedback.y,PI. u)
          annotation (Line(points={{19,136},{40,136}}, color={0,0,127}));
        connect(PI.y,limiter. u)
          annotation (Line(points={{63,136},{78.8,136}}, color={0,0,127}));
        connect(limiter.y, boundary_rear1.u0_var) annotation (Line(points={{92.6,136},{114,136},{114,182},{106,182}},
                                                                                                                    color={0,0,127}));
        connect(multiSensor_Tpm2.n_flow_out, feedback.u2) annotation (Line(points={{20,188},{10,188},{10,144}},
                                                                                                          color={0,0,127}));
        connect(multiSensor_Tpm1.fore, boundary_fore.rear) annotation (Line(
            points={{40,208},{86,208}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear.fore, flowResistanceA.rear) annotation (Line(
            points={{-74,208},{-68,208}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm.rear, flowResistanceA.fore) annotation (Line(
            points={{-40,208},{-48,208}},
            color={158,66,200},
            thickness=0.5));
        connect(ramp1.y, feedback.u1) annotation (Line(points={{-13,136},{2,136}}, color={0,0,127}));
        connect(boundary_fore1.u0_var, ramp3.y) annotation (Line(points={{-100,150},{-112,150},{-112,144},{-123,144}}, color={0,0,127}));
        connect(sinkA.inlet,multiSensor_Tpm5. outlet) annotation (Line(
            points={{-44,-64},{-32,-64}},
            color={158,66,200},
            thickness=0.5));
        connect(feedback2.y, PI2.u) annotation (Line(points={{-13,16},{-20,16}}, color={0,0,127}));
        connect(PI3.u,feedback3. y)
          annotation (Line(points={{-108,-64},{-115,-64}},
                                                         color={0,0,127}));
        connect(feedback3.u1,airFlow_setPoint2. y)
          annotation (Line(points={{-132,-64},{-138,-64},{-138,-48},{-139,-48}},
                                                       color={0,0,127}));
        connect(PI3.y,limiter3. u) annotation (Line(points={{-85,-64},{-77.2,-64}},
                                   color={0,0,127}));
        connect(sourceA.outlet, flowResistanceA1.inlet) annotation (Line(
            points={{102,-64},{90,-64}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm4.inlet, flowResistanceA1.outlet)
          annotation (Line(
            points={{60,-64},{70,-64}},
            color={158,66,200},
            thickness=0.5));
        connect(sinkB.inlet, flowResistanceB1.outlet) annotation (Line(
            points={{64,-8},{64,-16}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistanceB1.inlet, multiSensor_Tpm6.outlet)
          annotation (Line(
            points={{64,-36},{64,-48},{46,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(sourceB.outlet,multiSensor_Tpm7. inlet) annotation (Line(
            points={{-58,-30},{-58,-48},{-38,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(PI2.y, limiter2.u) annotation (Line(points={{-43,16},{-58,16},{-58,5.2}}, color={0,0,127}));
        connect(multiSensor_Tpm4.outlet, evaporator.inletA) annotation (Line(
            points={{38,-64},{20,-64},{20,-64},{12,-64}},
            color={158,66,200},
            thickness=0.5));
        connect(evaporator.outletA, multiSensor_Tpm5.inlet)
          annotation (Line(
            points={{-8,-64},{-12,-64}},
            color={158,66,200},
            thickness=0.5));
        connect(sinkA.u0_var,limiter3. y)
          annotation (Line(points={{-56,-64},{-63.4,-64}},
                                                         color={0,0,127}));
        connect(multiSensor_Tpm5.n_flow_out,feedback3. u2) annotation (Line(points={{-32,-68},{-42,-68},{-42,-82},{-124,-82},{-124,-72}},
                                                                 color={0,0,127}));
        connect(evaporator.outletB, multiSensor_Tpm6.inlet) annotation (Line(
            points={{12,-48},{26,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm7.outlet, evaporator.inletB) annotation (Line(
            points={{-18,-48},{-8,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm7.n_flow_out, feedback2.u2) annotation (Line(points={{-18,-44},{-4,-44},{-4,8}}, color={0,0,127}));
        connect(feedback2.u1, refFlow_setPoint1.y) annotation (Line(points={{4,16},{15,16}}, color={0,0,127}));
        connect(sourceB.u0_var, limiter2.y) annotation (Line(points={{-52,-18},{-52,-14},{-58,-14},{-58,-8.6}},
                                                                                            color={0,0,127}));
        connect(singleSensorSelect.inlet, evaporator.outletB)
          annotation (Line(
            points={{26,-24},{22,-24},{22,-48},{12,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(singleSensorSelect1.inlet, evaporator.inletB)
          annotation (Line(
            points={{-18,-24},{-14,-24},{-14,-48},{-8,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(sensorVaporQuality.inlet, evaporator.inletB)
          annotation (Line(
            points={{-18,-56},{-10,-56},{-10,-48},{-8,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(sensorVaporQuality1.inlet, evaporator.outletB)
          annotation (Line(
            points={{26,-56},{16,-56},{16,-48},{12,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(sinkA1.inlet, multiSensor_Tpm9.outlet) annotation (Line(
            points={{102,-166},{62,-166}},
            color={158,66,200},
            thickness=0.5));
        connect(feedback4.y, PI4.u) annotation (Line(points={{29,-244},{32,-244}}, color={0,0,127}));
        connect(PI5.u,feedback5. y)
          annotation (Line(points={{106,-130},{87,-130}},color={0,0,127}));
        connect(feedback5.u1,airFlow_setPoint3. y)
          annotation (Line(points={{70,-130},{29,-130}},
                                                       color={0,0,127}));
        connect(feedback5.u2,multiSensor_Tpm9. n_flow_out)
          annotation (Line(points={{78,-138},{78,-162},{62,-162}},color={0,0,127}));
        connect(sinkA1.u0_var, limiter5.y) annotation (Line(points={{114,-166},{127.4,-166}}, color={0,0,127}));
        connect(PI5.y,limiter5. u) annotation (Line(points={{129,-130},{146,-130},{146,-166},{141.2,-166}},
                                   color={0,0,127}));
        connect(sourceA1.outlet, flowResistanceA2.inlet) annotation (Line(
            points={{-102,-166},{-84,-166}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm8.inlet, flowResistanceA2.outlet)
          annotation (Line(
            points={{-48,-166},{-64,-166}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm11.outlet, condenser.inletB) annotation (Line(
            points={{20,-182},{8,-182}},
            color={158,66,200},
            thickness=0.5));
        connect(condenser.outletB, multiSensor_Tpm10.inlet)
          annotation (Line(
            points={{-12,-182},{-24,-182}},
            color={158,66,200},
            thickness=0.5));
        connect(condenser.inletA, multiSensor_Tpm8.outlet)
          annotation (Line(
            points={{-12,-166},{-14,-166},{-14,-166},{-26,-166}},
            color={158,66,200},
            thickness=0.5));
        connect(condenser.outletA, multiSensor_Tpm9.inlet) annotation (Line(
            points={{8,-166},{20,-166},{20,-166},{42,-166}},
            color={158,66,200},
            thickness=0.5));
        connect(sinkB1.inlet, flowResistanceB2.outlet) annotation (Line(
            points={{-76,-218},{-76,-212}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistanceB2.inlet, multiSensor_Tpm10.outlet)
          annotation (Line(
            points={{-76,-192},{-76,-182},{-44,-182}},
            color={158,66,200},
            thickness=0.5));
        connect(sourceB1.outlet, multiSensor_Tpm11.inlet)
          annotation (Line(
            points={{80,-196},{80,-182},{40,-182}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm11.n_flow_out, feedback4.u2) annotation (Line(points={{20,-186},{6,-186},{6,-236},{20,-236}}, color={0,0,127}));
        connect(PI4.y, limiter4.u) annotation (Line(points={{55,-244},{60.8,-244}}, color={0,0,127}));
        connect(refFlow_setPoint2.y, feedback4.u1) annotation (Line(points={{-3,-244},{12,-244}}, color={0,0,127}));
        connect(singleSensorSelect2.inlet, condenser.outletB)
          annotation (Line(
            points={{-24,-208},{-20,-208},{-20,-182},{-12,-182}},
            color={158,66,200},
            thickness=0.5));
        connect(singleSensorSelect3.inlet, condenser.inletB)
          annotation (Line(
            points={{22,-206},{14,-206},{14,-182},{8,-182}},
            color={158,66,200},
            thickness=0.5));
        connect(sensorVaporQuality2.inlet, condenser.outletB)
          annotation (Line(
            points={{-24,-220},{-20,-220},{-20,-182},{-12,-182}},
            color={158,66,200},
            thickness=0.5));
        connect(sensorVaporQuality3.inlet, condenser.inletB)
          annotation (Line(
            points={{22,-220},{14,-220},{14,-182},{8,-182}},
            color={158,66,200},
            thickness=0.5));
        connect(limiter4.y, sourceB1.u0_var) annotation (Line(points={{74.6,-244},{86,-244},{86,-208}}, color={0,0,127}));
        connect(rampChemicalPotential.y, boundary_rear.u0_var) annotation (Line(points={{-113,196},{-96,196},{-96,214},{-86,214}}, color={0,0,127}));
        connect(feedback1.u1, rampMassflow.y) annotation (Line(points={{42,240},{28,240},{28,244},{15,244}}, color={0,0,127}));
        connect(discretizedHEXUndir.rearA, multiSensor_Tpm.fore) annotation (Line(
            points={{-12,208},{-20,208}},
            color={158,66,200},
            thickness=0.5));
        connect(discretizedHEXUndir.foreB, multiSensor_Tpm3.rear) annotation (Line(
            points={{-12,192},{-22,192}},
            color={158,66,200},
            thickness=0.5));
        connect(discretizedHEXUndir.rearB, multiSensor_Tpm2.fore) annotation (Line(
            points={{8,192},{20,192}},
            color={158,66,200},
            thickness=0.5));
        connect(discretizedHEXUndir.foreA, multiSensor_Tpm1.rear) annotation (Line(
            points={{8,208},{20,208}},
            color={158,66,200},
            thickness=0.5));
        annotation (
          experiment(StopTime=30, Tolerance=1e-6, Interval=0.03, __Dymola_Algorithm="Dassl"),
          Icon(coordinateSystem(preserveAspectRatio=false)),
          Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-280,-300},{260,300}}),
            graphics={
              Rectangle(extent={{-208,-104},{170,-280}}, lineColor={158,66,200}),
              Text(
                extent={{-192,-106},{-106,-128}},
                textColor={158,66,200},
                textString="Condenser"),
              Text(
                extent={{-192,80},{-106,58}},
                textColor={158,66,200},
                textString="Evaporator"),
              Rectangle(extent={{-208,82},{170,-94}}, lineColor={158,66,200}),
              Text(
                extent={{-192,274},{-106,252}},
                textColor={158,66,200},
                textString="Undirected"),
              Rectangle(extent={{-208,276},{170,100}}, lineColor={158,66,200})}),
          Documentation(info="<html>
<u>Owner: <a href=\"mailto:niels.weber@dlr.de\">Niels Weber</a></u>
</html>"));
      end TestDiscretizedHEXvsDir;
    end Tests;

    package Internal
      extends Modelica.Icons.InternalPackage;

      model ConductionElementHEX "ConductionElement for single-phase fluids"
        extends PartialConductionElementHEX;

        parameter Modelica.Units.SI.CoefficientOfHeatTransfer U_nom=3000 "Heat transfer coefficient to medium";
        parameter Modelica.Units.SI.MolarFlowRate n_flow_nom=1 "Nominal mass-flow rate for heat transfer calculation"
          annotation (Dialog(group="Heat transfer parameters"));

        constant Real Re_exp = 0.8 "Reynolds-exponent for heat transfer calculation";

      initial equation
        assert(init <> Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.port,
          "This initialization will lead to large nonlinear equation systems. Please choose 'T0', 'h0', 'rear' or 'fore'.");

      equation
        //Estimation of heat transfer coefficient
        U = max(U_min, U_nom*(abs(n_flow/(n_flow_nom/nCellsParallel)))^Re_exp);

        annotation (Dialog(tab="Initialization", group="Enthalpy"),
                    Icon(coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<u>Undirected implementation of the Conduction Element for the DiscritizedHex.</u>
<u>Concerning the heat transfer it is assumed, that the main term influencing the coefficient of heat transfer is the mass flow rate. Therefore a nominal value for the heat transfer coefficient at a nominal mass flow rate can be set. Furthermore a minimum value U_min for the coefficient of heat transfer is set to ensure heat transfer at zero mass flow.</u>
<u>For further documentation see the documentation of the motherclass.</u>
</html>"));
      end ConductionElementHEX;

      model ConductionElementHEX_twoPhase "ConductionElement for two-phase fluids"
        extends PartialConductionElementHEX(redeclare replaceable package Medium =
              Chemical.Media.myMedia.Interfaces.PartialTwoPhaseMedium);

        import Modelica.Math;

        parameter Modelica.Units.SI.CoefficientOfHeatTransfer U_liq_nom=700 "Nominal coefficient of heat transfer for liquid region";
        parameter Modelica.Units.SI.CoefficientOfHeatTransfer U_vau_nom=500 "Nominal coefficient of heat transfer for vapour region";
        parameter Modelica.Units.SI.CoefficientOfHeatTransfer U_tu_nom=1000 "Nominal coefficient of heat transfer for two-phase region";

        parameter Modelica.Units.SI.MolarFlowRate n_flow_nom=0.3 "Nominal mass-flow rate for heat transfer calculation";


        constant Real Re_exu_cond(unit="1") = 0.4 "Reynolds-Exponent for heat transfer calculation at condensation (Yan&Lin, 1999)";
        constant Real Re_exu_evap(unit="1") = 0.5 "Reynolds-Exponent for heat transfer calculation at evaporation (Yan&Lin, 1999)";

        Real x "Vapor quality calculated from enthalpies";

      protected
        Modelica.Units.SI.CoefficientOfHeatTransfer U_liq "Coefficient of heat transfer for liquid region";
        Modelica.Units.SI.CoefficientOfHeatTransfer U_tp "Coefficient of heat transfer for two-phase region";
        Modelica.Units.SI.CoefficientOfHeatTransfer U_vap "Coefficient of heat transfer for vapour region";

        Modelica.Units.SI.MolarEnthalpy h_dew=Medium.dewEnthalpy(Medium.setSat_p(state)) "Dew enthalpy at inlet";
        Modelica.Units.SI.MolarEnthalpy h_bubble=Medium.bubbleEnthalpy(Medium.setSat_p(state)) "Bubble enthalpy at inlet";

      initial equation
        assert(init <> Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.port,
          "This initialization will lead to large nonlinear equation systems. Please choose 'h0', 'rear' or 'fore'.");
        assert(init <> Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.T,
          "Temperature might is not independent of pressrue in the 2-phase area. Please choose 'h0', 'rear' or 'fore'.");

      equation
        //1) Calculated from enthalpies --> can go below zero and above one (needed for heat transfer calculations)!
        x = (h - h_bubble)/(h_dew - h_bubble);

        U_liq = max(U_min, U_liq_nom*(abs(n_flow)/(n_flow_nom/nCellsParallel))^Re_exu_cond);
        U_vap = max(U_min, U_vau_nom*(abs(n_flow)/(n_flow_nom/nCellsParallel))^Re_exu_evap);
        U_tp = max(U_min, U_tu_nom);

        //Coefficient of heat transfer dependent on vapor quality (interpolation in phase-transition regions)
        U = smooth(1, noEvent(
            if x < -delta_x then U_liq
            elseif x < delta_x then U_liq + 0.5*(U_tp - U_liq)*(1 + Math.sin(x*Modelica.Constants.pi/(2*delta_x)))
            elseif x < 1 - delta_x then U_tp
            elseif x < 1 + delta_x then U_tp + 0.5*(U_vap - U_tp)*(1 + Math.sin((x - 1)*Modelica.Constants.pi/(2*delta_x)))
            else U_vap));

        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
             Line(
               points={{0,100},{0,-30}},
               color={238,46,47}),
             Line(
               points={{-100,0},{100,0}},
               thickness=0.5,
               color={158,66,200}),
             Ellipse(
               extent={{-70,-70},{70,70}},
               lineColor={158,66,200},
               lineThickness=0.5,
               fillColor={5,188,158},
               fillPattern=FillPattern.Solid,
               pattern=LinePattern.Solid,
                startAngle=45,
                endAngle=225),
             Line(
               points={{-50,-30},{50,-30}},
               color={238,46,47}),
             Line(
               points={{-50,-15},{50,-15}},
               color={238,46,47}),
             Line(
               points={{-50,0},{50,0}},
               color={238,46,47}),
             Line(
               points={{-50,15},{50,15}},
               color={238,46,47}),
             Line(
               points={{-50,30},{50,30}},
               color={238,46,47}),
             Ellipse(
               extent={{70,70},{-70,-70}},
               lineColor={158,66,200},
               lineThickness=0.5,
               pattern=LinePattern.Solid,
                startAngle=45,
                endAngle=225),
              Ellipse(
                extent={{-16,-54},{-6,-64}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{6,-36},{16,-46}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-42,-44},{-32,-54}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-18,-20},{-8,-30}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{6,6},{16,-4}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{54,0},{64,-10}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{38,34},{48,24}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{24,16},{34,6}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{20,-22},{30,-32}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{46,-16},{56,-26}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{34,2},{44,-8}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{2,-10},{12,-20}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{54,24},{64,14}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{14,-54},{24,-64}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{32,-38},{42,-48}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-26,-36},{-16,-46}},
                lineColor={158,66,200},
                fillColor={5,188,158},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<u>Undirected implementation of the Conduction Element for the DiscritizedHex.</u>
<u>Concerning the heat transfer coefficient it is assumed, that the main term influencing the coefficient of heat transfer is the mass flow rate. Therefore a nominal value for the heat transfer coefficient at a nominal mass flow rate can be set. The reynolds exponents for normalization of the heat transfer coefficient for evaporation and condensation are taken from Yan, Yi-Yie, &amp; Lin, T.-F. (1999). Condensation heat transfer and potential drop of refrigerant R-134a in a small pipe. International Journal of Heat and Mass Transfer, 42(4) and Yan, Y.-Y., &amp; Lin, T.-F. (1999). Evaporation Heat Transfer and ChemicalPotential Drop of Refrigerant R-134a in a Plate Heat Exchanger. Journal of Heat Transfer, 121(1). Furthermore a minimum value U_min for the coefficient of heat transfer is set to ensure heat transfer at zero mass flow.</u>
<u>For further documentation see the documentation of the motherclass.</u>
</html>"));
      end ConductionElementHEX_twoPhase;

      partial model PartialConductionElementHEX "Parent for CEs for discretizedHEX"
        extends Chemical.Undirected.Processes.Internal.PartialConductionElement(
                                                            final neglectChemicalPotentialChanges=true);

          parameter Modelica.Units.SI.Area A=1 "Contact area of volume with medium";

          parameter Integer nCellsParallel = 1 "Number of parallel discretization elements";

          constant Modelica.Units.SI.CoefficientOfHeatTransfer U_min=1 "Minimum heat transfer coefficient for temperature adaption at zero massflow";

          Modelica.Units.SI.CoefficientOfHeatTransfer U "Heat transfer coefficient to medium";
      equation
        k = U*A;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end PartialConductionElementHEX;

      partial model PartialDiscretizedHEX "Base class for undirected discretized heat exchangers"
        extends Chemical.HeatExchangers.Internal.DiscretizedHexIcon;

        replaceable package MediumA =
            Chemical.Media.myMedia.Interfaces.PartialMedium                           "Medium model side A" annotation (choicesAllMatching=true, Dialog(group="Medium definitions"));
        replaceable package MediumB =
            Chemical.Media.myMedia.Interfaces.PartialMedium                           "Medium model side B" annotation (choicesAllMatching=true, Dialog(group="Medium definitions"));

        replaceable model ConductionElementA =
            Chemical.Undirected.HeatExchangers.Internal.ConductionElementHEX                                    constrainedby
          Chemical.Undirected.HeatExchangers.Internal.PartialConductionElementHEX(
          final nCellsParallel=nCellsParallel,
          final A=A/nCells,
          final V=V_Hex/nCells,
          redeclare package Medium = MediumA,
          final enforce_global_energy_conservation=enforce_global_energy_conservation,
          final init=init_A,
          final h_0=h0_A) "Heat transfer element model for side A" annotation (choicesAllMatching=true, Dialog(group="Medium definitions"));
        replaceable model ConductionElementB =
            Chemical.Undirected.HeatExchangers.Internal.ConductionElementHEX                                    constrainedby
          Chemical.Undirected.HeatExchangers.Internal.PartialConductionElementHEX(
          final nCellsParallel=1,
          final A=A/nCells,
          final V=V_Hex/nCells,
          redeclare package Medium = MediumB,
          final enforce_global_energy_conservation=enforce_global_energy_conservation,
          final init=init_B,
          final h_0=h0_B) "Heat transfer element model for side B" annotation (choicesAllMatching=true, Dialog(group="Medium definitions"));

        parameter Boolean initializeMassFlow=false "Initialize mass flow at inlets?" annotation (Dialog(tab="Initialization", group="Mass flow"));
        parameter Modelica.Units.SI.MolarFlowRate n_flow_0_A=0 "Initial mass flow for side A"
          annotation (Dialog(
            tab="Initialization",
            group="Mass flow",
            enable=initializeMassFlow));
        parameter Modelica.Units.SI.MolarFlowRate n_flow_0_B=0 "Initial mass flow for side B"
          annotation (Dialog(
            tab="Initialization",
            group="Mass flow",
            enable=initializeMassFlow));
        parameter Integer nCells=3 "Number of discretization elements";
        parameter Modelica.Units.SI.Area A=10 "Conductive area of heat exchanger" annotation (Dialog(group="Heat transfer parameters"));
        parameter Modelica.Units.SI.Volume V_Hex=0.001 "Volume for heat transfer calculation" annotation (Dialog(group="Heat transfer parameters"));
        parameter Boolean enforce_global_energy_conservation=false "If true, exact global energy conservation is enforced by feeding back all energy stored locally back in the system" annotation (Dialog(tab="Advanced"));
        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg
          "Regularization mass flow to switch between positive- and negative-massflow model" annotation (Dialog(tab="Advanced"));

        parameter Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement init_A=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.h "Initialization method for h side A" annotation (Dialog(tab="Initialization", group="Enthalpy"), choices(
            choice=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.h "h0",
            choice=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.rear "rear",
            choice=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.fore "fore"));
        parameter Modelica.Units.SI.MolarEnthalpy h0_A=MediumA.h_default "Initial enthalpy side A" annotation (Dialog(
            tab="Initialization",
            group="Enthalpy",
            enable=(init_A == Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.h)));
        parameter Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement init_B=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.h "Initialization method for h side B" annotation (Dialog(tab="Initialization", group="Enthalpy"), choices(
            choice=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.h "h0",
            choice=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.rear "rear",
            choice=Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.fore "fore"));
        parameter Modelica.Units.SI.MolarEnthalpy h0_B=MediumB.h_default "Initial enthalpy side B" annotation (Dialog(
            tab="Initialization",
            group="Enthalpy",
            enable=(init_B == Chemical.Undirected.Processes.Internal.InitializationMethodsCondElement.h)));

        //Parameterization of HEX Wall
        parameter Modelica.Units.SI.CoefficientOfHeatTransfer k_wall=100 "Coefficient of heat transfer for pipe wall" annotation (Dialog(group="Heat transfer parameters"));

        parameter Boolean calculate_efficiency=false "Enable calculation of efficiency";

        Modelica.Units.SI.HeatFlowRate Q_flow_A=sum(thermalElementA.heatPort.Q_flow);
        Modelica.Units.SI.HeatFlowRate Q_flow_B=sum(thermalElementB.heatPort.Q_flow);
        Modelica.Units.SI.MolarFlowRate n_flow_A=rearA.n_flow;
        Modelica.Units.SI.MolarFlowRate n_flow_B=rearB.n_flow;
        Modelica.Units.SI.Mass M_A=sum(thermalElementA.M);
        Modelica.Units.SI.Mass M_B=sum(thermalElementB.M);
        Modelica.Units.SI.Energy deltaE_system=sum(thermalElementA.deltaE_system) + sum(thermalElementB.deltaE_system);

        Chemical.HeatExchangers.Internal.DiscretizedHEXSummary summary "Summary record of quantities";

      protected
        parameter Boolean crossFlow=false "Selection whether HEX is in crossflow or counterflow configuration";
        parameter Integer nCellsParallel=1 "Number of discretization elements in parallel";
        parameter Modelica.Units.SI.ThermalConductance G=k_wall*A "Wall thermal conductance" annotation (Dialog(group="Wall parameters"));
        outer Chemical.DropOfCommons dropOfCommons;

        function efficiency =
            Chemical.HeatExchangers.Internal.calculateEfficiency (                  redeclare
              package                                                                                 MediumA = MediumA, redeclare
              package                                                                                                                      MediumB = MediumB);

        // no regstep since this is only used as a output
        MediumA.ThermodynamicState stateA_in=if noEvent(rearA.n_flow) > 0 then rearA.state_forwards else foreA.state_rearwards;
        MediumA.ThermodynamicState stateA_out=if noEvent(rearA.n_flow) > 0 then foreA.state_forwards else rearA.state_rearwards;
        MediumB.ThermodynamicState stateB_in=if noEvent(rearB.n_flow) > 0 then rearB.state_forwards else foreB.state_rearwards;
        MediumB.ThermodynamicState stateB_out=if noEvent(rearB.n_flow) > 0 then foreB.state_forwards else rearB.state_rearwards;

      public
        Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor[nCells](each G=G/nCells) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={0,0})));

        ConductionElementA thermalElementA[nCells] annotation (Placement(transformation(extent={{10,-90},{-10,-70}})));
        ConductionElementB thermalElementB[nCells] annotation (Placement(transformation(extent={{-10,90},{10,70}})));

        Chemical.Undirected.Interfaces.Rear rearA(redeclare package Medium = MediumA) annotation (Placement(transformation(extent={{110,-90},{90,-70}}),
              iconTransformation(extent=if crossFlow then {{-110,-10},{-90,10}} else {{110,-90},{90,-70}}, rotation=if crossFlow then -90 else 0)));
        Chemical.Undirected.Interfaces.Fore foreA(redeclare package Medium = MediumA) annotation (Placement(transformation(extent={{-90,-90},{-110,-70}}),
              iconTransformation(extent=if crossFlow then {{90,-10},{110,10}} else {{-90,-90},{-110,-70}}, rotation=if crossFlow then -90 else 0)));
        Chemical.Undirected.Interfaces.Rear rearB(redeclare package Medium = MediumB) annotation (Placement(transformation(extent={{-110,70},{-90,90}}),
              iconTransformation(extent=if crossFlow then {{110,-90},{90,-70}} else {{-110,70},{-90,90}})));
        Chemical.Undirected.Interfaces.Fore foreB(redeclare package Medium = MediumB) annotation (Placement(transformation(extent={{90,70},{110,90}}),
              iconTransformation(extent=if crossFlow then {{-90,-90},{-110,-70}} else {{90,70},{110,90}})));

      equation
        //Summary record
        summary.Tin_A = MediumA.temperature(stateA_in);
        summary.Tin_B = MediumB.temperature(stateB_in);
        summary.Tout_A = MediumA.temperature(stateA_out);
        summary.Tout_B = MediumB.temperature(stateB_out);
        summary.hin_A = MediumA.specificEnthalpy(stateA_in);
        summary.hin_B = MediumB.specificEnthalpy(stateB_in);
        summary.hout_A = MediumA.specificEnthalpy(stateA_out);
        summary.hout_B = MediumB.specificEnthalpy(stateB_out);
        summary.dT_A = summary.Tout_A - summary.Tin_A;
        summary.dT_B = summary.Tout_B - summary.Tin_B;
        summary.dh_A = summary.hout_A - summary.hin_A;
        summary.dh_B = summary.hout_B - summary.hin_B;

        summary.efficiency = if calculate_efficiency then efficiency(
          stateA_in,
          stateB_in,
          stateA_out,
          stateB_out,
          rearA.n_flow,
          rearB.n_flow,
          Q_flow_A) else 0;

        summary.Q_flow_A = Q_flow_A;
        summary.Q_flow_B = Q_flow_B;

        annotation (Documentation(info="<html>
<u>
This is the partial parent class for undirected discretized heat exchangers. It contains
the rear and fore connectors as well as a number of conduction elements (which
is set by the parameter <code>nCells</code>) as discrete control volumes to
exchange heat between two fluid streams.
</u>
<u>
The conduction elements are computing a heat transfer coefficient between their
heatport and the fluid contained. They are replaceable with a choice between
a single-phase and a two-phase version, both can be further parametrized.
Although the single-phase version works for two-phase media (not the other
way around), using the two-phase one for two-phase media enables to set different
heat transfer coefficients depending on the phase (liquid/gaseous/2-phase) state
of the medium.
</u>
<u>
Note that since the model uses conductionElements as discrete control volumes that
in turn assume quasi-stationary mass and, therefore, are part of a fluid stream
rather than break it into two (like a full volume would), the same holds for
both sides of the heat exchanger &ndash; they are part of a fluid stream and
don&apos;t break it. The quasi-stationary mass assumption also implies that for
(fast) changing masses/densities in any of the conduction elements the heat
exchanger will (slightly) violate the conservation of energy. Furthermore, the
conduction elements change their behavior for reversed mass flow, therefore, this
model asserts for negative mass flow with the level
&quot;<a href=\"Chemical.DropOfCommons\">DropOfCommons</a>.assertionLevel&quot;.
</u>
<u>
The parameters <code>A</code> (heat transferring area), <code>k_wall</code> (heat
transfer coefficient of the wall between the streams) and the heat transfer
coefficients in the conduction elements scale the transferred heat (the middle
only if the wall and the latter only of the heat transfer into a fluid is the
choke of the heatflow).
</u>
<u>
The parameter <code>V</code> determines the amount of fluid in the heat exchanger
and, therefore, the dynamic for non-steady states.
</u>
<u>
The &quot;Initialization&quot; tab allows for a mass flow initialization for both
paths, as well as to determine from which direction the enthalpy in the control
volumes should be initialized (fore/rear), or if it should start with a given
enthalpy. The other option is to initialize the enthalpy with a given value.
</u>
<u>
The &quot;Advanced&quot; tab allows to modify the mass flow that triggers the
reverse-mass-flow-assertion and has an option to enforce global conservation of
energy. The latter is done by feeding back any energy the conduction elements
accumulated over time, basically making it impossible to store energy in their
fluid long-term. While this enforces long-term conservation of energy, it
changes the medium-/short-term dynamics of the system and is, therefore,
disabled by default.
</u>
</html>"));
      end PartialDiscretizedHEX;
    end Internal;
    annotation(Icon(graphics={
        Line(
          points={{-66,58}},
          color={158,66,200},
          thickness=0.5),
          Polygon(
            points={{76,60},{-84,-40},{-78,-22},{-98,-22},{12,46},{0,60},{76,60}},
            pattern=LinePattern.None,
            lineThickness=1,
            fillColor={158,66,200},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-76,-60},{84,40},{78,22},{98,22},{-12,-46},{0,-60},{-76,-60}},
            pattern=LinePattern.None,
            lineThickness=1,
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid,
            lineColor={0,0,0})}),
                        Documentation(info="<html>
<u>This package contains several types of heat exchangers with different model approaches.</u>
</html>",   revisions="<html>
<u><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</u>
</html>"));
  end HeatExchangers;

  package FlowControl "Package for undirected flow control components"
    extends Modelica.Icons.Package;

    model BasicControlValve
      "Basic valve model with optional flow characteristics for incompressible fluids"
      extends Chemical.Undirected.FlowControl.Internal.PartialValve;

      import FlowCoeffType =
             Chemical.FlowControl.Internal.Types.FlowCoefficientTypesBasic;

      replaceable function valveCharacteristics =
        Chemical.FlowControl.Internal.ControlValve.linearCharacteristics
        constrainedby Chemical.FlowControl.Internal.ControlValve.partialValveCharacteristics
        "Select valve characteristics"
        annotation (
          choicesAllMatching = true,
          Dialog(group = "Valve parameters"),
          Documentation(info="<html>
<u>Characteristic curve of the valve.</u>
</html>"));

      parameter FlowCoeffType flowCoefficient = FlowCoeffType.Kvs "Select type of flow coefficient" annotation(Dialog(group = "Valve parameters"));
      //Reference Values
      parameter Real Kvs(unit = "m3/h")  "Kvs-value (metric) from data sheet (valve fully open)" annotation(Evaluate = true,
        Dialog(group = "Valve parameters",enable = (flowCoefficient ==FlowCoeffType.
              Kvs)));
      parameter Real Cvs_US "Cvs-value (US [gal/min]) from data sheet (valve fully open)" annotation(Evaluate = true,
      Dialog(group = "Valve parameters",enable = (flowCoefficient ==FlowCoeffType.Cvs_US)));
      parameter Real Cvs_UK "Cvs-value (UK [gal/min]) from data sheet (valve fully open)" annotation(Evaluate = true,
      Dialog(group = "Valve parameters",enable = (flowCoefficient ==FlowCoeffType.Cvs_UK)));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_ref_set "Set reference mass flow in kg/s"
        annotation (Evaluate=true, Dialog(group="Valve parameters", enable=(flowCoefficient == FlowCoeffType.n_flow_set)));

    protected
      Modelica.Units.SI.VolumeFlowRate V_flow_ref=if flowCoefficient == FlowCoeffType.Kvs then Kvs/secondsPerHour elseif flowCoefficient == FlowCoeffType.Cvs_US
           then (Cvs_US/1.1561)/secondsPerHour elseif flowCoefficient == FlowCoeffType.Cvs_UK then (Cvs_UK/0.9626)/secondsPerHour else n_flow_ref_set/rho_ref
        "Reference volume flow";

    equation
      //Calculate reference mass flow from reference volume flow
      n_flow_ref = V_flow_ref*rho_ref;

      k_u = valveCharacteristics(u, k_min);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-84,0},{-40,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{40,0},{-40,0}},
              color={158,66,200},
              thickness=0.5,
              pattern=LinePattern.Dash),
            Line(
              points={{0,0},{0,60}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{40,0},{84,0}},
              color={158,66,200},
              thickness=0.5),
            Polygon(
              points={{-20,40},{0,0},{20,40},{-20,40}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor=DynamicSelect({255,255,255}, if invertInput == true then
                      {158,66,200} else {255,255,255}),
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-20,20},{0,-20},{20,20},{-20,20}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor=DynamicSelect({255,255,255}, if invertInput == true then
                      {158,66,200} else {255,255,255}),
              fillPattern=FillPattern.Solid,
              origin={0,-20},
              rotation=180)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Undirected implementation of the Basic Control Valve.</u>
<u>This model serves for most incompressible applications where basic control valves are needed. </u>
<u><br>The modeler has the ability to choose between different valve characteristics and flow coefficients.</u>
<u>The three standard curve characteristics (linear, parabolic, equal-percentage) are implemented and can be chosen.</u>
<u><br>To conclude the parameterization, a flow coefficient has to be set. Most data sheets of valves deliver a corresponding &quot;KVs (CVs)&quot;-Value. Otherwise a nominal mass-flow rate can be set. </u>
<u>For incompressible flow, the reference values for density (1g/cm3) and potential (1bar) should be unchanged.</u>
</html>"));
    end BasicControlValve;

    model SpecificValveType "Specific technical valve types"
      extends Chemical.Undirected.FlowControl.Internal.PartialValve;

      import FlowCoeffType =
             Chemical.FlowControl.Internal.Types.FlowCoefficientTypes;

      replaceable record ZetaValueRecord =
          Chemical.FlowControl.Internal.Curves.SlideValveZetaCurve
        constrainedby
        Chemical.FlowControl.Internal.Curves.PartialCharacteristicZetaCurves "Select valve type"
          annotation(choicesAllMatching = true, Dialog(group = "Valve parameters"));

      parameter FlowCoeffType flowCoefficient = FlowCoeffType.Kvs "Select type of flow coefficient" annotation(Dialog(group = "Valve parameters"));
      //Set valve data as parameter
      parameter Modelica.Units.SI.Diameter d_valve "Flow diameter" annotation (
          Evaluate=true, Dialog(group="Valve parameters", enable=(flowCoefficient
               == FlowCoeffType.flowDiameter)));
      //Reference Values
      parameter Real Kvs(unit = "m3/h")  "Kvs-value (metric) from data sheet (valve fully open)" annotation(Evaluate = true,
        Dialog(group = "Valve parameters",enable = (flowCoefficient ==FlowCoeffType.
              Kvs)));
      parameter Real Cvs_US "Cvs-value (US [gal/min]) from data sheet (valve fully open)" annotation(Evaluate = true,
      Dialog(group = "Valve parameters",enable = (flowCoefficient ==FlowCoeffType.Cvs_US)));
      parameter Real Cvs_UK "Cvs-value (UK [gal/min]) from data sheet (valve fully open)" annotation(Evaluate = true,
      Dialog(group = "Valve parameters",enable = (flowCoefficient ==FlowCoeffType.Cvs_UK)));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_ref_set "Set reference mass flow in kg/s"
        annotation (Evaluate=true, Dialog(group="Valve parameters", enable=(flowCoefficient == FlowCoeffType.n_flow_set)));

    protected
      constant ZetaValueRecord valveData;
      Modelica.Units.SI.Area A_valve=0.25*Modelica.Constants.pi*d_valve^2
        "Cross-sectional valve area";

      Real k_u(unit="1") "Kv/Kvs, respecting flow characteristics";
      Real k_u_zeta(unit="1") "Kv/Kvs respecting zeta curve";

      Modelica.Blocks.Tables.CombiTable1Ds combiTable1D_zeta(
      final smoothness = Modelica.Blocks.Types.Smoothness.MonotoneContinuousDerivative1,
      final tableOnFile = false,
      final table = valveData.zetaTable) "Interpolation of zeta datapoints";

      Real zeta(unit="1", start = 0) "zeta value for potential loss calculation";
      Real zeta1(unit="1") = valveData.zetaTable[end,2] "zeta value for fully open valve";

      Modelica.Units.SI.VolumeFlowRate V_flow_ref=if flowCoefficient == FlowCoeffType.Kvs then Kvs/secondsPerHour elseif flowCoefficient == FlowCoeffType.Cvs_US
           then (Cvs_US/1.1561)/secondsPerHour elseif flowCoefficient == FlowCoeffType.Cvs_UK then (Cvs_UK/0.9626)/secondsPerHour elseif flowCoefficient ==
          FlowCoeffType.flowDiameter then A_valve*sqrt((2/zeta1)*(du_ref/rho_ref)) else n_flow_ref_set/rho_ref "Reference volume flow";

    equation
      //Calculate reference mass flow from reference volume flow
      n_flow_ref = V_flow_ref*rho_ref;

      //Retrieving zeta value from actuation signal
      combiTable1D_zeta.u = u;
      zeta = combiTable1D_zeta.y[1];

      //Evaluate characteristic for given zeta curve
      k_u_zeta = sqrt(zeta1/zeta);

      k_u = k_min + (1 - k_min)*k_u_zeta;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-84,0},{-40,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{40,0},{-40,0}},
              color={158,66,200},
              thickness=0.5,
              pattern=LinePattern.Dash),
            Line(
              points={{0,0},{0,60}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{40,0},{84,0}},
              color={158,66,200},
              thickness=0.5),
            Polygon(
              points={{-20,40},{0,0},{20,40},{-20,40}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor=DynamicSelect({255,255,255}, if invertInput == true then
                      {158,66,200} else {255,255,255}),
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-20,20},{0,-20},{20,20},{-20,20}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor=DynamicSelect({255,255,255}, if invertInput == true then
                      {158,66,200} else {255,255,255}),
              fillPattern=FillPattern.Solid,
              origin={0,-20},
              rotation=180)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Undirected implementation of the specific valve type.</u>
<u>This valve models the behavior of specific valve types.</u>
<u><br>The technical type of the valve can be chosen (e.g. sliding valve). The characteristic curve is then set accordingly from a table for the zeta (flow resistance) values dependent on the valve opening.</u>
<u><br>To conclude the parameterization, a flow coefficient has to be set. Most data sheets of valves deliver a corresponding &quot;KVs (CVs)&quot;-Value. Otherwise a nominal mass-flow rate or a flow-diameter can be set. </u>
<u>For incompressible flow, the reference values for density (1g/cm3) and potential (1bar) should be unchanged.</u>
</html>"));
    end SpecificValveType;

    model TanValve "Valve with tan-shaped flow resistance"
      extends Chemical.Undirected.Interfaces.SISOBiFlow(final cliu_u_out=
            true);

      Modelica.Blocks.Interfaces.RealInput u(unit="1") "Valve control input []"
        annotation (Placement(
            transformation(extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,80}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,80})));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the flow" annotation (Dialog(tab="Advanced"));
      parameter Boolean invertInput = false "Zero represents a closed valve for non-inverted, open for inverted";
      parameter Modelica.Units.SI.MolarFlowRate n_flow_ref=0.1 "Reference mass flow";
      parameter Modelica.Units.SI.ChemicalPotential u_ref=1e5 "Reference potential";
      parameter Real relativeLeakiness(unit="1") = 1e-3 "Imperfection of valve";

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      Real k(unit="(Pa.s)/kg");
      Real u2(unit="1");

    equation
      if invertInput then
        u2 = max(relativeLeakiness,min(1-relativeLeakiness,u));
      else
        u2 = max(relativeLeakiness,min(1-relativeLeakiness,1-u));
      end if;

      k = u_ref/n_flow_ref*tan(u2*Modelica.Constants.pi/2);

      //forwards model
      du_fore = -k*n_flow;
      h_fore_out = h_rear_in;
      Xi_fore_out = Xi_rear_in;

      //rearwards model
      du_rear = - du_fore; // - because of the inverted direction
      h_rear_out = h_fore_in;
      Xi_rear_out = Xi_fore_in;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-84,0},{-40,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{40,0},{-40,0}},
              color={158,66,200},
              thickness=0.5,
              pattern=LinePattern.Dash),
            Line(
              points={{0,0},{0,60}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{40,0},{84,0}},
              color={158,66,200},
              thickness=0.5),
            Polygon(
              points={{-20,40},{0,0},{20,40},{-20,40}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor=DynamicSelect({255,255,255}, if invertInput == true then
                      {158,66,200} else {255,255,255}),
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-20,20},{0,-20},{20,20},{-20,20}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor=DynamicSelect({255,255,255}, if invertInput == true then
                      {158,66,200} else {255,255,255}),
              fillPattern=FillPattern.Solid,
              origin={0,-20},
              rotation=180)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Undirected implementation of the Tan valve. </u>
<u>The TanValve is the most basic valve and can be used when no valve type is set yet. </u>
<u>It adjusts its flow resistance coefficient according to a tangens of the input. The pole of the tan function can lead to numerical problems.</u>
</html>"));
    end TanValve;

    model MCV "Massflow control valve"
      extends Chemical.Undirected.Interfaces.SISOBiFlow(
        final cliu_u_out=false,
        final L=100,
        final u_min=u_min_par);

      import Mode = Chemical.FlowControl.Internal.Types.MassflowControlValveMode;

      Modelica.Blocks.Interfaces.RealInput setpoint_var if setpointFromInput "Desired mass-flow [kg/s or m3/s]"
        annotation (Placement(
            transformation(extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,80}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,80})));

      Modelica.Blocks.Interfaces.RealOutput clippingOutput = (dp - du_int) if enableClippingOutput ""
        annotation (Placement(
            transformation(extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-80}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-80})));

      parameter Mode mode = Mode.mass_flow "Valve mode";
      parameter Boolean setpointFromInput = false "= true, if desired flow rate is set via setpoint_var input";
      parameter Modelica.Units.SI.MolarFlowRate nassFlow_set_par=0 "Mass flow variable to set"
        annotation (Dialog(enable=(not setpointFromInput) and mode == Mode.mass_flow));
      parameter Modelica.Units.SI.VolumeFlowRate volumeFlow_set_par=0 "Mass flow variable to set"
        annotation (Dialog(enable=(not setpointFromInput) and mode == Mode.volume_flow));
      parameter Modelica.Units.SI.Time TC=0.1 "Time constant of setpoint dynamic";
      parameter Real k1(unit="1") = 100 "Timeconstant factor"
        annotation(Dialog(tab="Advanced"));
      parameter Real k2(unit="1") = 100 "Integrator windup factor"
        annotation(Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.ChemicalPotential u_min_par=dropOfCommons.u_min "Minimal steady-state output potential" annotation (Dialog(tab="Advanced"));
      parameter Boolean enableClippingOutput = false "= true, if clippingOutput enabled";

      Modelica.Units.SI.Density rho_rear_in=Medium.density(rear.state_forwards);
      Modelica.Units.SI.Density rho_fore_in=Medium.density(fore.state_rearwards);
      Modelica.Units.SI.Density rho_in=Chemical.Undirected.Internal.regStep(
              n_flow,
              rho_rear_in,
              rho_fore_in,
              n_flow_reg);

      Modelica.Units.SI.VolumeFlowRate V_flow_fore=n_flow/rho_rear_in;
      Modelica.Units.SI.VolumeFlowRate V_flow_rear=n_flow/rho_fore_in;
      Modelica.Units.SI.VolumeFlowRate V_flow=Chemical.Undirected.Internal.regStep(
              n_flow,
              V_flow_fore,
              V_flow_rear,
              n_flow_reg);

      constant Modelica.Units.SI.ChemicalPotential eps=1;
      Modelica.Units.SI.ChemicalPotential dp=Chemical.Undirected.Internal.regStep(
              n_flow,
              du_fore,
              -du_rear,
              n_flow_reg);

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      Modelica.Units.SI.MolarFlowRate n_flow_set;
      Modelica.Units.SI.VolumeFlowRate V_flow_set;
      Modelica.Blocks.Interfaces.RealInput setpoint "Internal setpoint connector";

      Modelica.Units.SI.ChemicalPotential dr=fore.r - rear.r;
      Modelica.Units.SI.ChemicalPotential dr_set;
    public
      Modelica.Units.SI.ChemicalPotential du_int(start=-1e5);
      Modelica.Units.SI.ChemicalPotential du_corr_fore=k2*(du_fore - du_int);
      Modelica.Units.SI.ChemicalPotential du_corr_rear=k2*(du_rear + du_int);
      Modelica.Units.SI.ChemicalPotential du_corr=Chemical.Undirected.Internal.regStep(
              n_flow,
              du_corr_fore,
              -du_corr_rear,
              n_flow_reg);

    initial equation
      du_int = 0;

    equation
      connect(setpoint_var, setpoint);
      if setpointFromInput then
        n_flow_set = setpoint;
        V_flow_set = setpoint;
      else
        setpoint = 0;
        n_flow_set = massFlow_set_par;
        V_flow_set = volumeFlow_set_par;
      end if;

      // compute dr_set required for desired mass-flow or Volume-flow dynamic
      if mode==Mode.mass_flow then
        dr_set = - L/TC * (n_flow_set - n_flow);
      else
        dr_set = - L/TC*(V_flow_set*rho_in - n_flow);
      end if;

      // compute potential drop dynamic very fast, so dr tracks dr_set.
      // dr is limited, since it can be very high for non-smooth systems (e.g. a jump in input potential)
      TC/k1 * der(du_int) = max(-1e8, min(1e8,dr)) - dr_set + du_corr;

      // limit dp to a so that u_out > u_min and no potential is created
      du_fore = max(u_min - u_rear_in, min(0, du_int));
      du_rear = max(u_min - u_fore_in, min(0, -du_int));

      h_fore_out = h_rear_in;
      h_rear_out = h_fore_in;

      Xi_fore_out = Xi_rear_in;
      Xi_rear_out = Xi_fore_in;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Ellipse(
              extent={{-56,54},{64,-66}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{-84,0},{84,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{40,0},{-48,0}},
              color={158,66,200},
              thickness=0.5,
              pattern=LinePattern.Dash),
            Line(
              points={{-52,-30},{52,-30}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{-52,30},{52,30}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,60}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{40,-60},{60,-80}},
              lineColor={0,0,0},
              fillColor = DynamicSelect({255,255,255}, if abs(dp - du_int) <= eps then {0,140,72} else {238,46,47}),
              fillPattern=FillPattern.Solid)}),
                  Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Undirected implementation of MCV.</u>
<u>This component can be used to emulate a mass-flow regulated valve. </u>
<u>It works similar to the directed MCV.</u>
</html>"));
    end MCV;

    model CheckValve "Valve that allows only positive mass_flow"
      extends Chemical.Undirected.Interfaces.SISOBiFlow(final cliu_u_out=
            false);

      parameter Modelica.Units.SI.MolarFlowRate n_flow_ref=dropOfCommons.n_flow_reg "Reference mass flow" annotation (Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.ChemicalPotential u_ref=1e5 "Reference potential" annotation (Dialog(tab="Advanced"));

    equation
      //forwards model
      du_fore = if n_flow < 0 then u_ref*((n_flow/n_flow_ref)^2) else 0;
      h_fore_out = h_rear_in;
      Xi_fore_out = Xi_rear_in;

      //rearwards model
      du_rear = - du_fore; // - because of the inverted direction
      h_rear_out = h_fore_in;
      Xi_rear_out = Xi_fore_in;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Ellipse(
              extent={{-56,54},{64,-66}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{-84,0},{84,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{40,0},{-48,0}},
              color={158,66,200},
              thickness=0.5,
              pattern=LinePattern.Dash),
            Line(
              points={{-52,-30},{52,-30}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{-52,30},{52,30}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,30},{20,10}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,-30},{20,-10}},
              color={158,66,200},
              thickness=0.5)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Undirected implementation of CheckValve.</u>
<u>Valve that allows positive mass_flow and builds up a large potential difference against negative mass_flow.</u>
</html>"));
    end CheckValve;

    package Tests "Tests for undirected FlowControl components"
    extends Modelica.Icons.ExamplesPackage;

      model CheckValve "Test for undirected CheckValve"
        extends Modelica.Icons.Example;

        replaceable package Medium = Chemical.Media.myMedia.Air.SimpleAir constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium package"
          annotation (choicesAllMatching=true, Documentation(info="<html>
<u>
Medium package used in the Test.
</u>
</html>"));

        inner Chemical.DropOfCommons dropOfCommons(assertionLevel=AssertionLevel.warning)
          annotation (Placement(transformation(extent={{52,-82},{72,-62}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300) annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=200000)
          annotation (Placement(transformation(extent={{52,-10},{72,10}})));
        FlowControl.CheckValve checkValve(
          redeclare package stateOfMatter = stateOfMatter,
          L=1e-4,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state) annotation (Placement(transformation(extent={{-22,-10},{-2,10}})));
        Modelica.Blocks.Sources.Pulse pulse(
          amplitude=2e5,
          period=0.5,
          offset=1e5)
          annotation (Placement(transformation(extent={{-96,-10},{-76,10}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.1,
          l=10,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss)
          annotation (Placement(transformation(extent={{14,-10},{34,10}})));
      equation
        connect(pulse.y, boundary_rear.u0_var) annotation (Line(points={{-75,0},{-60,0},{-60,6},{-42,6}}, color={0,0,127}));

        connect(boundary_rear.fore, checkValve.rear) annotation (Line(
            points={{-30,0},{-22,0}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance.rear, checkValve.fore) annotation (Line(
            points={{14,0},{-2,0}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore.rear, flowResistance.fore) annotation (Line(
            points={{52,0},{34,0}},
            color={158,66,200},
            thickness=0.5));
        annotation (
          experiment(StopTime=1, Tolerance=1e-6, Interval=0.001, __Dymola_Algorithm="Dassl"),
          Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
              Documentation(info="<html>
<u>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
      end CheckValve;

      model TanValve "Test for undirected TanValve"
        extends Modelica.Icons.Example;

        replaceable package Medium = Chemical.Media.myMedia.Air.SimpleAir constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium package"
          annotation (choicesAllMatching=true, Documentation(info="<html>
<u>
Medium package used in the Test.
</u>
</html>"));

        inner Chemical.DropOfCommons dropOfCommons(assertionLevel=AssertionLevel.warning)
          annotation (Placement(transformation(extent={{52,-82},{72,-62}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300) annotation (Placement(transformation(extent={{-106,-10},{-86,10}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=100000)
          annotation (Placement(transformation(extent={{96,-10},{116,10}})));
        FlowControl.TanValve tanValve(redeclare package stateOfMatter =
              stateOfMatter,                                                           invertInput=false) annotation (Placement(transformation(extent={{16,26},{36,46}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-52,-10},{-32,10}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore1(redeclare
            package                                                                  stateOfMatter =
              stateOfMatter,                                                                                        u0_par=100000)
          annotation (Placement(transformation(extent={{96,26},{116,46}})));
        Chemical.Undirected.Topology.JunctionRFF2 junctionRFF2_1(redeclare
            package                                                                stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-24,-10},{-4,10}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm(redeclare
            package                                                                   stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-78,-2},{-58,18}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm1(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{26,-2},{46,18}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm2(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-12,34},{8,54}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm3(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{42,34},{62,54}})));
        Modelica.Blocks.Sources.Ramp ramp(
          height=1,
          duration=2,
          offset=0,
          startTime=50)
          annotation (Placement(transformation(extent={{-12,64},{8,84}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance1(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{68,26},{88,46}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance2(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{64,-10},{84,10}})));
        Modelica.Blocks.Sources.Pulse pulse(
          amplitude=4e4,
          period=10,
          offset=0.8e5)
          annotation (Placement(transformation(extent={{-138,-10},{-118,10}})));
      equation
        connect(ramp.y, tanValve.u)
          annotation (Line(points={{9,74},{26,74},{26,44}}, color={0,0,127}));
        connect(flowResistance1.rear, multiSensor_Tpm3.fore) annotation (Line(
            points={{68,36},{62,36}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore1.rear, flowResistance1.fore) annotation (Line(
            points={{96,36},{88,36}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore.rear, flowResistance2.fore) annotation (Line(
            points={{96,0},{84,0}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance2.rear, multiSensor_Tpm1.fore) annotation (Line(
            points={{64,0},{46,0}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm1.rear, junctionRFF2_1.foreB) annotation (Line(
            points={{26,0},{-4,0}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm2.fore, tanValve.rear) annotation (Line(
            points={{8,36},{16,36}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm2.rear, junctionRFF2_1.foreA)
          annotation (Line(
            points={{-12,36},{-14,36},{-14,10}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionRFF2_1.rear, flowResistance.fore) annotation (Line(
            points={{-24,0},{-32,0}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance.rear, multiSensor_Tpm.fore) annotation (Line(
            points={{-52,0},{-58,0}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm.rear, boundary_rear.fore) annotation (Line(
            points={{-78,0},{-86,0}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm3.rear, tanValve.fore) annotation (Line(
            points={{42,36},{36,36}},
            color={158,66,200},
            thickness=0.5));
        connect(pulse.y, boundary_rear.u0_var) annotation (Line(points={{-117,0},{-108,0},{-108,6},{-98,6}},
                                                                                           color={0,0,127}));
        annotation (Diagram(coordinateSystem(extent={{-120,-100},{120,100}})),
          experiment(StopTime=100, Tolerance=1e-6, Interval=0.1, __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
        <u>Owner: <a href=\"mailto:niels.weber@dlr.de\">Niels Weber</a></u>
</html>"));
      end TanValve;

      model BasicControlValve "Test for undirected BasicControlValve"
        extends Modelica.Icons.Example;

        replaceable package Medium =
            Chemical.Media.myMedia.Water.ConstantPropertyLiquidWater
          constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium package"
          annotation (choicesAllMatching=true, Documentation(info="<html>
<u>
Medium package used in the Test.
</u>
</html>"));

        inner Chemical.DropOfCommons dropOfCommons(assertionLevel = AssertionLevel.warning)
          annotation (Placement(transformation(extent={{-170,94},{-150,114}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300) annotation (Placement(transformation(extent={{-116,50},{-96,70}})));
        FlowControl.BasicControlValve valveLinear(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          flowCoefficient=Chemical.FlowControl.Internal.Types.FlowCoefficientTypesBasic.Kvs,
          Kvs=5,
          redeclare function valveCharacteristics =
              Chemical.FlowControl.Internal.ControlValve.linearCharacteristics)
          annotation (Placement(transformation(extent={{-10,50},{10,70}})));

        Chemical.Undirected.Processes.FlowResistance flowResistance(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{-80,50},{-60,70}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore1(redeclare
            package                                                                  stateOfMatter =
              stateOfMatter,                                                                                        u0_par=100000)
          annotation (Placement(transformation(extent={{96,50},{116,70}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm2(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-40,58},{-20,78}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm3(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{26,58},{46,78}})));
        Modelica.Blocks.Sources.Ramp ramp(
          height=1,
          duration=10,
          offset=0,
          startTime=5)
          annotation (Placement(transformation(extent={{174,10},{154,30}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear2(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300) annotation (Placement(transformation(extent={{-116,-10},{-96,10}})));
        FlowControl.BasicControlValve valveParabolic(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          invertInput=false,
          flowCoefficient=Chemical.FlowControl.Internal.Types.FlowCoefficientTypesBasic.Kvs,
          Kvs=5,
          redeclare function valveCharacteristics =
              Chemical.FlowControl.Internal.ControlValve.parabolicCharacteristics,
          n_flow_ref_set=0.1) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore2(redeclare
            package                                                                  stateOfMatter =
              stateOfMatter,                                                                                        u0_par=100000)
          annotation (Placement(transformation(extent={{96,-10},{116,10}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm4(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-40,-2},{-20,18}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm5(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{26,-2},{46,18}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear1(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300) annotation (Placement(transformation(extent={{-116,-70},{-96,-50}})));
        FlowControl.BasicControlValve valveEqualPercentage(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          flowCoefficient=Chemical.FlowControl.Internal.Types.FlowCoefficientTypesBasic.Kvs,
          Kvs=5,
          redeclare function valveCharacteristics =
              Chemical.FlowControl.Internal.ControlValve.equalPercentageCharacteristics)
          annotation (Placement(transformation(extent={{-10,-70},{10,-50}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=100000)
          annotation (Placement(transformation(extent={{96,-70},{116,-50}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm7(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-40,-62},{-20,-42}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm8(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{26,-62},{46,-42}})));
        Modelica.Blocks.Sources.Pulse pulse(
          amplitude=2e4,
          period=2,
          offset=0.9e5)
          annotation (Placement(transformation(extent={{-172,-10},{-152,10}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance6(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{60,50},{80,70}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance1(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{60,-10},{80,10}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance2(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance3(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{-80,-70},{-60,-50}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance4(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{60,-70},{80,-50}})));
      equation
        connect(boundary_rear1.u0_var, pulse.y) annotation (Line(points={{-108,-54},{-140,-54},{-140,0},{-151,0}}, color={0,0,127}));
        connect(ramp.y, valveParabolic.u_in)
          annotation (Line(points={{153,20},{0,20},{0,8}}, color={0,0,127}));
        connect(ramp.y, valveEqualPercentage.u_in) annotation (Line(points={{153,20},
                {120,20},{120,-40},{0,-40},{0,-52}}, color={0,0,127}));
        connect(boundary_rear.u0_var, pulse.y) annotation (Line(points={{-108,66},{-140,66},{-140,0},{-151,0}}, color={0,0,127}));
        connect(pulse.y, boundary_rear2.u0_var) annotation (Line(points={{-151,0},{-130,0},{-130,6},{-108,6}},
                                                                                             color={0,0,127}));
        connect(valveLinear.u_in, ramp.y) annotation (Line(points={{0,68},{0,90},{120,
                90},{120,20},{153,20}}, color={0,0,127}));
        connect(flowResistance.rear, boundary_rear.fore) annotation (Line(
            points={{-80,60},{-96,60}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance.fore, multiSensor_Tpm2.rear) annotation (Line(
            points={{-60,60},{-40,60}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm2.fore, valveLinear.rear) annotation (Line(
            points={{-20,60},{-10,60}},
            color={158,66,200},
            thickness=0.5));
        connect(valveLinear.fore, multiSensor_Tpm3.rear)
          annotation (Line(
            points={{10,60},{18,60},{18,60},{26,60}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm3.fore, flowResistance6.rear) annotation (Line(
            points={{46,60},{60,60}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance6.fore, boundary_fore1.rear) annotation (Line(
            points={{80,60},{96,60}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore2.rear, flowResistance1.fore) annotation (Line(
            points={{96,0},{80,0}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance1.rear, multiSensor_Tpm5.fore) annotation (Line(
            points={{60,0},{46,0}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm5.rear, valveParabolic.fore) annotation (Line(
            points={{26,0},{10,0}},
            color={158,66,200},
            thickness=0.5));
        connect(valveParabolic.rear, multiSensor_Tpm4.fore) annotation (Line(
            points={{-10,0},{-20,0}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm4.rear, flowResistance2.fore) annotation (Line(
            points={{-40,0},{-60,0}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance2.rear, boundary_rear2.fore) annotation (Line(
            points={{-80,0},{-96,0}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear1.fore, flowResistance3.rear) annotation (Line(
            points={{-96,-60},{-80,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance3.fore, multiSensor_Tpm7.rear)
          annotation (Line(
            points={{-60,-60},{-40,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm7.fore, valveEqualPercentage.rear)
          annotation (Line(
            points={{-20,-60},{-10,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(valveEqualPercentage.fore, multiSensor_Tpm8.rear)
          annotation (Line(
            points={{10,-60},{26,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm8.fore, flowResistance4.rear) annotation (Line(
            points={{46,-60},{60,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance4.fore, boundary_fore.rear) annotation (Line(
            points={{80,-60},{96,-60}},
            color={158,66,200},
            thickness=0.5));
        annotation (Diagram(coordinateSystem(extent={{-180,-100},{180,120}})),
          experiment(
            StopTime=20,
         Tolerance=1e-6,
         Interval=0.02,
            __Dymola_Algorithm="Dassl"),
          Icon(coordinateSystem),
          Documentation(info="<html>
        <u>Owner: <a href=\"mailto:niels.weber@dlr.de\">Niels Weber</a></u>
</html>"));
      end BasicControlValve;

      model SpecificValveType "Test for undirected SpecificValveType"
        extends Modelica.Icons.Example;

        replaceable package Medium =
            Chemical.Media.myMedia.Water.ConstantPropertyLiquidWater
          constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium package"
          annotation (choicesAllMatching=true, Documentation(info="<html>
<u>
Medium package used in the Test.
</u>
</html>"));

        inner Chemical.DropOfCommons dropOfCommons(assertionLevel = AssertionLevel.warning)
          annotation (Placement(transformation(extent={{-170,94},{-150,114}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear1(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300) annotation (Placement(transformation(extent={{-116,50},{-96,70}})));
        FlowControl.SpecificValveType slideValve(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          redeclare record ZetaValueRecord =
              Chemical.FlowControl.Internal.Curves.SlideValveZetaCurve,
          flowCoefficient=Chemical.FlowControl.Internal.Types.FlowCoefficientTypes.Kvs,
          Kvs=5) annotation (Placement(transformation(extent={{-10,50},{10,70}})));

        Chemical.Undirected.Processes.FlowResistance flowResistance(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{-80,50},{-60,70}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore1(redeclare
            package                                                                  stateOfMatter =
              stateOfMatter,                                                                                        u0_par=100000)
          annotation (Placement(transformation(extent={{96,50},{116,70}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm2(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-40,58},{-20,78}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm3(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{26,58},{46,78}})));
        Modelica.Blocks.Sources.Ramp ramp(
          height=1,
          duration=10,
          offset=0,
          startTime=5)
          annotation (Placement(transformation(extent={{174,10},{154,30}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300) annotation (Placement(transformation(extent={{-116,-10},{-96,10}})));
        FlowControl.SpecificValveType slideValveInverse(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          invertInput=true,
          redeclare record ZetaValueRecord =
              Chemical.FlowControl.Internal.Curves.SlideValveZetaCurve,
          flowCoefficient=Chemical.FlowControl.Internal.Types.FlowCoefficientTypes.Kvs,
          d_valve=0.005,
          Kvs=5,
          n_flow_ref_set=0.1) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=100000)
          annotation (Placement(transformation(extent={{96,-10},{116,10}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm4(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-40,-2},{-20,18}})));
        Chemical.Undirected.Sensors.MultiSensor_Tpm multiSensor_Tpm5(redeclare
            package                                                                    stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{26,-2},{46,18}})));

        Chemical.Undirected.Processes.FlowResistance flowResistance6(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{60,50},{80,70}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance1(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{60,-10},{80,10}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance2(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.05,
          l=1,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
              (                                                                                                       k=1e3))
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Modelica.Blocks.Sources.Pulse pulse(
          amplitude=2e4,
          period=2,
          offset=0.9e5)
          annotation (Placement(transformation(extent={{-170,-10},{-150,10}})));
      equation
        connect(ramp.y, slideValveInverse.u_in)
          annotation (Line(points={{153,20},{0,20},{0,8}}, color={0,0,127}));
        connect(slideValve.u_in, ramp.y) annotation (Line(points={{0,68},{0,90},{120,
                90},{120,20},{153,20}}, color={0,0,127}));
        connect(pulse.y, boundary_rear.u0_var) annotation (Line(points={{-149,0},{-128,0},{-128,6},{-108,6}},
                                                                                            color={0,0,127}));
        connect(boundary_rear1.u0_var, boundary_rear.u0_var) annotation (Line(points={{-108,66},{-130,66},{-130,6},{-108,6}}, color={0,0,127}));
        connect(boundary_rear1.fore, flowResistance.rear)
          annotation (Line(
            points={{-96,60},{-90,60},{-90,60},{-80,60}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance2.rear, boundary_rear.fore) annotation (Line(
            points={{-80,0},{-96,0}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance.fore, multiSensor_Tpm2.rear) annotation (Line(
            points={{-60,60},{-40,60}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm4.rear, flowResistance2.fore) annotation (Line(
            points={{-40,0},{-60,0}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm4.fore, slideValveInverse.rear) annotation (Line(
            points={{-20,0},{-10,0}},
            color={158,66,200},
            thickness=0.5));
        connect(slideValveInverse.fore, multiSensor_Tpm5.rear)
          annotation (Line(
            points={{10,0},{26,0},{26,0}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm5.fore, flowResistance1.rear) annotation (Line(
            points={{46,0},{60,0}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance1.fore, boundary_fore.rear) annotation (Line(
            points={{80,0},{96,0}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore1.rear, flowResistance6.fore) annotation (Line(
            points={{96,60},{80,60}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance6.rear, multiSensor_Tpm3.fore) annotation (Line(
            points={{60,60},{46,60}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm3.rear, slideValve.fore) annotation (Line(
            points={{26,60},{10,60}},
            color={158,66,200},
            thickness=0.5));
        connect(slideValve.rear, multiSensor_Tpm2.fore)
          annotation (Line(
            points={{-10,60},{-20,60},{-20,60}},
            color={158,66,200},
            thickness=0.5));
        annotation (Diagram(coordinateSystem(extent={{-180,-100},{180,120}})),
          experiment(
            StopTime=20,
         Tolerance=1e-6,
         Interval=0.02,
            __Dymola_Algorithm="Dassl"),
          Icon(coordinateSystem),
          Documentation(info="<html>
        <u>Owner: <a href=\"mailto:niels.weber@dlr.de\">Niels Weber</a></u>
</html>"));
      end SpecificValveType;

      model MCV "Test for undirected MCV"
        extends Modelica.Icons.Example;

        replaceable package Medium = Chemical.Media.myMedia.Air.SimpleAir constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium package"
            annotation (choicesAllMatching=true, Documentation(info="<html>
<u>
Medium package used in the Test.
</u>
</html>"));

        inner Chemical.DropOfCommons dropOfCommons(assertionLevel=AssertionLevel.warning)
          annotation (Placement(transformation(extent={{100,0},{120,20}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear2(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=false,
          T0_par(displayUnit="K") = 300,
          u0_par=200000) annotation (Placement(transformation(extent={{-30,30},{-10,50}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore(redeclare
            package                                                                stateOfMatter =
              stateOfMatter,                                                                                      u0_par=100000)
          annotation (Placement(transformation(extent={{60,30},{80,50}})));
        Chemical.Undirected.FlowControl.MCV mCV(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
          setpointFromInput=true,
          massFlow_set_par=0.1,
          volumeFlow_set_par=1) annotation (Placement(transformation(extent={{0,30},{20,50}})));

        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear3(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300,
          u0_par=200000) annotation (Placement(transformation(extent={{-30,0},{-10,20}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore2(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=100000)
          annotation (Placement(transformation(extent={{60,0},{80,20}})));
        Chemical.Undirected.FlowControl.MCV mCV1(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
          massFlow_set_par=0.1,
          volumeFlow_set_par=1) annotation (Placement(transformation(extent={{0,0},{20,20}})));

        Modelica.Blocks.Sources.Pulse pulse1(
          amplitude=2e5,
          period=2.5,
          offset=0.5e5,
          startTime=0)
          annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore6(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300,
          u0_par=200000) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-20,-20})));

        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear1(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=100000)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={70,-20})));
        Chemical.Undirected.FlowControl.MCV mCV2(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
          massFlow_set_par=0.1,
          volumeFlow_set_par=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=180,
              origin={10,-20})));

        Modelica.Blocks.Sources.Trapezoid trapezoid(
          amplitude=2e5,
          rising=0.5,
          width=0.75,
          falling=0.5,
          period=2.5,
          offset=0.5e5,
          startTime=0)
          annotation (Placement(transformation(extent={{-70,-30},{-50,-10}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore8(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par(displayUnit="K") = 300,
          u0_par=200000) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-20,-50})));

        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear5(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=100000)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={70,-50})));
        Chemical.Undirected.FlowControl.MCV mCV3(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
          massFlow_set_par=0.1,
          volumeFlow_set_par=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=180,
              origin={10,-50})));

        Modelica.Blocks.Sources.Pulse pulse3(
          amplitude=2e5,
          period=2.5,
          offset=0.5e5,
          startTime=0)
          annotation (Placement(transformation(extent={{-70,-60},{-50,-40}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore7(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=false,
          T0_par(displayUnit="K") = 300,
          u0_par=100000) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-20,-80})));

        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear4(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par=573.15,
          u0_par=100000) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={70,-80})));
        Chemical.Undirected.FlowControl.MCV mCV4(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
          massFlow_set_par=0.1,
          volumeFlow_set_par=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=180,
              origin={10,-80})));

        Modelica.Blocks.Sources.Trapezoid
                                      trapezoid1(
          amplitude=2e5,
          rising=0.5,
          width=1.5,
          falling=0.5,
          period=2.5,
          offset=0.5e5,
          startTime=0)
          annotation (Placement(transformation(extent={{112,-88},{92,-68}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=false,
          T0_par(displayUnit="K") = 300,
          u0_par=200000) annotation (Placement(transformation(extent={{-30,70},{-10,90}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore1(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=100000)
          annotation (Placement(transformation(extent={{60,70},{80,90}})));
        Chemical.Undirected.FlowControl.MCV mCV5(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.volume_flow,
          setpointFromInput=true,
          massFlow_set_par=0.1,
          volumeFlow_set_par=1) annotation (Placement(transformation(extent={{0,70},{20,90}})));

        Modelica.Blocks.Sources.Pulse pulse5(
          amplitude=0.2,
          period=2.5,
          offset=0.1,
          startTime=0)
          annotation (Placement(transformation(extent={{-70,50},{-50,70}})));
        Modelica.Blocks.Sources.Trapezoid trapezoid2(
          amplitude=5,
          rising=0.5,
          width=0.75,
          falling=0.5,
          period=2.5,
          offset=1,
          startTime=0)
          annotation (Placement(transformation(extent={{-70,90},{-50,110}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance5(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.1,
          l=10,
          L_value=0.01,
          computeL=false,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss)
          annotation (Placement(transformation(extent={{28,70},{48,90}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.1,
          l=10,
          L_value=0.01,
          computeL=false,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss)
          annotation (Placement(transformation(extent={{30,30},{50,50}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance1(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.1,
          l=10,
          L_value=0.01,
          computeL=false,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss)
          annotation (Placement(transformation(extent={{28,0},{48,20}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance2(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.1,
          l=10,
          L_value=0.01,
          computeL=false,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={40,-20})));
        Chemical.Undirected.Processes.FlowResistance flowResistance3(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.1,
          l=10,
          L_value=0.01,
          computeL=false,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={40,-50})));
        Chemical.Undirected.Processes.FlowResistance flowResistance4(
          redeclare package stateOfMatter = stateOfMatter,
          r=0.1,
          l=10,
          L_value=0.01,
          computeL=false,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={40,-80})));
        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear6(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=false,
          T0_par(displayUnit="K") = 300,
          u0_par=200000) annotation (Placement(transformation(extent={{-32,158},{-12,178}})));
        Chemical.Undirected.FlowControl.MCV mCV6(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
          setpointFromInput=false,
          massFlow_set_par=1,
          volumeFlow_set_par=1) annotation (Placement(transformation(extent={{-2,158},{18,178}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore3(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=200000)
          annotation (Placement(transformation(extent={{30,158},{50,178}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear7(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=false,
          T0_par(displayUnit="K") = 300,
          u0_par=500000) annotation (Placement(transformation(extent={{-32,116},{-12,136}})));
        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore4(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=0.001)
          annotation (Placement(transformation(extent={{30,116},{50,136}})));
        Chemical.Undirected.FlowControl.MCV mCV7(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
          setpointFromInput=false,
          massFlow_set_par=1,
          volumeFlow_set_par=1) annotation (Placement(transformation(extent={{-2,116},{18,136}})));

        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore10(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=false,
          T0_par(displayUnit="K") = 300,
          u0_par=200000) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={64,168})));
        Chemical.Undirected.FlowControl.MCV mCV8(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
          setpointFromInput=false,
          massFlow_set_par=1,
          volumeFlow_set_par=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=180,
              origin={94,168})));

        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear8(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=200000)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={126,168})));
        Chemical.Undirected.Boundaries.BoundaryFore boundaryFore11(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=false,
          T0_par(displayUnit="K") = 300,
          u0_par=500000) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={64,126})));
        Chemical.Undirected.Boundaries.BoundaryRear boundaryRear9(redeclare
            package                                                                 stateOfMatter =
              stateOfMatter,                                                                                       u0_par=0.001)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={126,126})));
        Chemical.Undirected.FlowControl.MCV mCV9(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          mode=Chemical.FlowControl.Internal.Types.MassflowControlValveMode.mass_flow,
          setpointFromInput=false,
          massFlow_set_par=1,
          volumeFlow_set_par=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=180,
              origin={94,126})));

      equation

        connect(pulse1.y, boundaryRear3.u0_var) annotation (Line(points={{-49,10},{-30,
                10},{-30,16},{-22,16}}, color={0,0,127}));
        connect(trapezoid.y, boundaryFore6.u0_var) annotation (Line(points={{-49,-20},
                {-30,-20},{-30,-26},{-22,-26}}, color={0,0,127}));
        connect(pulse3.y, boundaryFore8.u0_var) annotation (Line(points={{-49,-50},{-30,
                -50},{-30,-56},{-22,-56}}, color={0,0,127}));
        connect(pulse5.y,mCV. setpoint_var) annotation (Line(points={{-49,60},{10,60},
                {10,48}}, color={0,0,127}));
        connect(trapezoid2.y,mCV5. setpoint_var) annotation (Line(points={{-49,100},{10,
                100},{10,88}}, color={0,0,127}));
        connect(boundaryRear.fore, mCV5.rear) annotation (Line(
            points={{-10,80},{0,80}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear2.fore, mCV.rear) annotation (Line(
            points={{-10,40},{0,40}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear3.fore, mCV1.rear) annotation (Line(
            points={{-10,10},{0,10}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryFore6.rear, mCV2.fore) annotation (Line(
            points={{-10,-20},{-1.77636e-15,-20}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryFore8.rear, mCV3.fore) annotation (Line(
            points={{-10,-50},{-1.77636e-15,-50}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryFore7.rear, mCV4.fore) annotation (Line(
            points={{-10,-80},{0,-80}},
            color={158,66,200},
            thickness=0.5));
        connect(mCV5.fore, flowResistance5.rear) annotation (Line(
            points={{20,80},{28,80}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance5.fore, boundaryFore1.rear) annotation (Line(
            points={{48,80},{60,80}},
            color={158,66,200},
            thickness=0.5));
        connect(mCV.fore, flowResistance.rear) annotation (Line(
            points={{20,40},{30,40}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance.fore, boundaryFore.rear) annotation (Line(
            points={{50,40},{60,40}},
            color={158,66,200},
            thickness=0.5));
        connect(mCV1.fore, flowResistance1.rear) annotation (Line(
            points={{20,10},{28,10}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance1.fore, boundaryFore2.rear) annotation (Line(
            points={{48,10},{60,10}},
            color={158,66,200},
            thickness=0.5));
        connect(mCV2.rear, flowResistance2.fore) annotation (Line(
            points={{20,-20},{30,-20}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance2.rear, boundaryRear1.fore) annotation (Line(
            points={{50,-20},{60,-20}},
            color={158,66,200},
            thickness=0.5));
        connect(mCV3.rear, flowResistance3.fore) annotation (Line(
            points={{20,-50},{30,-50}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance3.rear, boundaryRear5.fore) annotation (Line(
            points={{50,-50},{60,-50}},
            color={158,66,200},
            thickness=0.5));
        connect(mCV4.rear, flowResistance4.fore) annotation (Line(
            points={{20,-80},{30,-80}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance4.rear, boundaryRear4.fore) annotation (Line(
            points={{50,-80},{60,-80}},
            color={158,66,200},
            thickness=0.5));
        connect(trapezoid1.y, boundaryRear4.u0_var) annotation (Line(points={{91,-78},
                {84,-78},{84,-86},{72,-86}}, color={0,0,127}));
        connect(boundaryRear6.fore, mCV6.rear) annotation (Line(
            points={{-12,168},{-2,168}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear7.fore, mCV7.rear) annotation (Line(
            points={{-12,126},{-2,126}},
            color={158,66,200},
            thickness=0.5));
        connect(mCV6.fore, boundaryFore3.rear) annotation (Line(
            points={{18,168},{30,168}},
            color={158,66,200},
            thickness=0.5));
        connect(mCV7.fore, boundaryFore4.rear) annotation (Line(
            points={{18,126},{30,126}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryFore10.rear, mCV8.fore) annotation (Line(
            points={{74,168},{84,168}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryFore11.rear, mCV9.fore) annotation (Line(
            points={{74,126},{84,126}},
            color={158,66,200},
            thickness=0.5));
        connect(mCV8.rear, boundaryRear8.fore) annotation (Line(
            points={{104,168},{116,168}},
            color={158,66,200},
            thickness=0.5));
        connect(boundaryRear9.fore, mCV9.rear) annotation (Line(
            points={{116,126},{104,126}},
            color={158,66,200},
            thickness=0.5));
        annotation (
          experiment(StopTime=10, Tolerance=1e-6, Interval=0.001, __Dymola_Algorithm="Dassl"),
          Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
              Documentation(info="<html>
<u>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
      end MCV;
    annotation (Documentation(info="<html>
<u>Tests&nbsp;for&nbsp;undirected&nbsp;FlowControl&nbsp;components</u>
</html>"));
    end Tests;

    package Internal "Internal models, functions and reckords for undirected FlowControl"

      extends Modelica.Icons.InternalPackage;

      partial model PartialValve "Partial implementation of a physical valve"
        extends Chemical.Undirected.Interfaces.SISOBiFlow(
                                      final cliu_u_out=true);

        parameter Boolean invertInput = false "Non-inverted: 0=closed, 1=open";
        parameter Real k_min(unit="1", min = 0.001, max = 1) = 0.03 "Remaining flow at actuation signal u = 0 (fraction of maximum mass flow at u = 1)";

        Modelica.Blocks.Interfaces.RealInput u_in(unit="1") "Valve control signal []"
          annotation (Placement(
              transformation(
              extent={{-20,-20},{20,20}},
              rotation=270,
              origin={0,80})));

        Real u(unit="1") "actuation input for flow calculation";

        parameter Modelica.Units.SI.ChemicalPotential du_ref=1e5
          "Reference potential difference"
          annotation (Dialog(tab="Advanced", group="Reference values"));
        parameter Modelica.Units.SI.Density rho_ref=1000 "Reference density"
          annotation (Dialog(tab="Advanced", group="Reference values"));

      protected
        constant Real secondsPerHour(final unit="s/h") = 3600 "Parameter for unit conversion";

        //Medium properties
        Modelica.Units.SI.Density rho_rear_in=Medium.density(rear.state_forwards);
        Modelica.Units.SI.Density rho_fore_in=Medium.density(fore.state_rearwards);

        Modelica.Units.SI.MolarFlowRate n_flow_ref "Reference mass flow derived from flow coefficient inputs";
        Real k_u(unit="1") "Kv/Kvs, respecting flow characteristics";

      equation
        //Inversion of input signal, actuation has to be between 0 and 1
        assert(u_in <=1 and u_in >=0, "Actuator signal out of valid range [0,1]", level = AssertionLevel.warning);
        if invertInput then
            u = max(0, min(1, 1-u_in));
        else
            u = max(0, min(1, u_in));
        end if;

         //forwards model
        du_fore = -sign(n_flow)*du_ref*(rho_ref/rho_rear_in)*(1/(k_u)*(n_flow/n_flow_ref))^2;
        h_fore_out = h_rear_in;
        Xi_fore_out = Xi_rear_in;

        //rearwards model
        du_rear = sign(n_flow)*du_ref*(rho_ref/rho_fore_in)*(1/(k_u)*(n_flow/n_flow_ref))^2;
        h_rear_out = h_fore_in;
        Xi_rear_out = Xi_fore_in;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<u>Partial&nbsp;implementation&nbsp;of&nbsp;a&nbsp;physical&nbsp;valve.</u>
</html>"));
      end PartialValve;
    annotation (Documentation(info="<html>
<u>Internal&nbsp;models,&nbsp;functions&nbsp;and&nbsp;reckords&nbsp;for&nbsp;undirected&nbsp;FlowControl</u>
</html>"));
    end Internal;
  annotation (Documentation(revisions="<html>
<u><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</u>

</html>",   info="<html>
<u>Undirected implementation of the flow control components. </u><u>These include physical Valves, as well as fow-control valves, that are not physical models, but rather numerical boundaries on the flow.</u>
</html>"),   Icon(graphics={
          Line(
            points={{-94,0},{-40,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{40,0},{-40,0}},
            color={158,66,200},
            thickness=0.5,
            pattern=LinePattern.Dash),
          Line(
            points={{40,0},{94,0}},
            color={158,66,200},
            thickness=0.5),
          Polygon(
            points={{-20,40},{0,0},{20,40},{-20,40}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor=DynamicSelect({255,255,255}, if invertInput == true then
                    {158,66,200} else {255,255,255}),
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-20,20},{0,-20},{20,20},{-20,20}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor=DynamicSelect({255,255,255}, if invertInput == true then
                    {158,66,200} else {255,255,255}),
            fillPattern=FillPattern.Solid,
            origin={0,-20},
            rotation=180)}));
  end FlowControl;

  package Sensors "Sensors package for undirected chemical simulation"
    extends Modelica.Icons.SensorsPackage;

    import Chemical.Sensors.Internal;

    model SingleSensorSelect "Sensor with selectable measured quantity"
      extends Internal.PartialSensor;
      import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

      parameter Chemical.Sensors.Internal.Types.Quantities quantity "Quantity to be measured";
      parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimum density" annotation (Dialog(tab="Advanced", group="Regularization"));
      parameter Boolean outputValue = false "Enable sensor-value output"
        annotation(Dialog(group="Output Value"));
      parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
        annotation(Dialog(group="Output Value", enable=outputValue));
      parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
        annotation(Dialog(tab="Initialization", enable=filter_output));
      parameter Real value_0(unit=Chemical.Sensors.Internal.getUnit(quantity)) = 0 "Initial output state of sensor"
        annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
      parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));

      Modelica.Blocks.Interfaces.RealOutput value_out(unit=Chemical.Sensors.Internal.getUnit(quantity)) = value if outputValue "Measured quantity [variable]"
        annotation (Placement(
            transformation(extent={{80,-20},{120,20}}),
              iconTransformation(extent={{80,-20},{120,20}})));

      function getQuantity = Chemical.Sensors.Internal.getQuantity (
        redeclare package stateOfMatter = stateOfMatter) "Quantity compute function"
        annotation (Documentation(info="<html>
      <u>This function computes the selected quantity from state. r and rho_min are neddet for the quantities r/u_total and v respectively.</u>
      </html>"));

      Real value(unit=Chemical.Sensors.Internal.getUnit(quantity));

    protected
      Real direct_value(unit=Chemical.Sensors.Internal.getUnit(quantity));

    initial equation
      if filter_output and init==InitMode.steadyState then
        value= direct_value;
      elseif filter_output then
        value = value_0;
      end if;

    equation

      direct_value = getQuantity(state, rear.r, quantity, rho_min);

      if filter_output then
        der(value) * TC = direct_value-value;
      else
        value = direct_value;
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-54,24},{66,-36}},
              lineColor={0,0,0},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{0,0},{0,-80}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-5,-75},{5,-85}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-60,30},{60,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-60,30},{60,-30}},
              textColor={158,66,200},
              textString=DynamicSelect("value", String(value, format="1."+String(digits)+"f"))),
            Text(
              extent={{2,17},{62,67}},
              textColor={175,175,175},
              textString="%quantity")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Undirected sensor for a single selectable quantity. For some quatities several units are available.</u>
</html>"));
    end SingleSensorSelect;

    model TwoPhaseSensorSelect "Sensor for a selectable quantity of a twoPhaseMedium"
      extends Internal.PartialSensor(redeclare package Medium=Medium2Phase);

      import Quantities=Chemical.Sensors.Internal.Types.TwoPhaseQuantities;
      import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

      replaceable package Medium2Phase =
          Chemical.Media.myMedia.Interfaces.PartialTwoPhaseMedium
                                                         "Medium model"
        annotation (choicesAllMatching=true,
          Documentation(info="<html>
<u>Replaceable medium package for the sensor. Medium must be a TwoPase Medium.</u>
</html>"));

      parameter Quantities quantity "Quantity the sensor measures"
        annotation(choicesAllMatching=true);
      parameter Boolean outputValue = false "Enable sensor-value output"
        annotation(Dialog(group="Output Value"));
      parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
        annotation(Dialog(group="Output Value", enable=outputValue));
      parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
        annotation(choicesAllMatching=true, Dialog(tab="Initialization", enable=filter_output));
      parameter Real value_0(unit=Chemical.Sensors.Internal.getTwoPhaseUnit(quantity)) = 0 "Initial output state of sensor"
        annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
      parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));

      Modelica.Blocks.Interfaces.RealOutput value_out(unit=Chemical.Sensors.Internal.getTwoPhaseUnit(quantity)) = value if outputValue "Measured quantity [variable]"
        annotation (Placement(
            transformation(extent={{80,-20},{120,20}}),
              iconTransformation(extent={{80,-20},{120,20}})));

      Real value(unit=Chemical.Sensors.Internal.getTwoPhaseUnit(quantity));

    protected
      Real direct_value(unit=Chemical.Sensors.Internal.getTwoPhaseUnit(quantity));

      function getQuantity = Chemical.Sensors.Internal.getTwoPhaseQuantity (
        redeclare package stateOfMatter = stateOfMatter) "Quantity compute function"
        annotation (Documentation(info="<html>
    <u>This function computes the selected two-phase quantity from state.</u>
      </html>"));

    initial equation
      if filter_output and init==InitMode.steadyState then
        value= direct_value;
      elseif filter_output then
        value = value_0;
      end if;

    equation
      direct_value = getQuantity(state, quantity);

      if filter_output then
        der(value) * TC = direct_value-value;
      else
        value = direct_value;
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-54,24},{66,-36}},
              lineColor={0,0,0},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{-100,-80},{100,-80}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,-80}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-5,-75},{5,-85}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-60,30},{60,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-60,30},{60,-30}},
              textColor={158,66,200},
              textString=DynamicSelect("value", String(value, format="1."+String(digits)+"f"))),
            Text(
              extent={{0,19},{60,69}},
              textColor={175,175,175},
              textString="%quantity")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Undirected sensor for a vapor quality. It must be separate from SingleSensorSelect, because it needs a TwoPhaseMedium.</u>
</html>"));
    end TwoPhaseSensorSelect;

    model SensorState "Sensor for whole state"
      extends Internal.PartialSensor;

      replaceable package Medium =
          Chemical.Media.myMedia.Interfaces.PartialMedium
        "Medium model"
        annotation (choicesAllMatching=true,
          Documentation(info="<html>
        <u>Medium Model for the sensor. Make sure it is the same as for all lines the sensors input is connected.</u>
        </html>"));

      Chemical.Interfaces.StateOutput state_out(redeclare package stateOfMatter =
            stateOfMatter)                                                                       "Measured value [variable]"
        annotation (Placement(transformation(extent={{80,-20},{120,20}})));

    equation

      state_out.state = state;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-54,24},{66,-36}},
              lineColor={0,0,0},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{0,0},{0,-80}},
              color={158,66,200},
              thickness=0.5),
            Rectangle(
              extent={{-60,30},{60,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-60,30},{60,-30}},
              textColor={158,66,200},
              textString="state"),
            Ellipse(
              extent={{-5,-75},{5,-85}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Sensor for measuring the full state.</u>
<u>This sensor can be connected to a fluid stream without a junction.</u>
</html>"));
    end SensorState;

    model SingleSensorX "Sensor for mass fraction of mixture"
      extends Internal.PartialSensor;

      import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

      replaceable package Medium =
          Chemical.Media.myMedia.Interfaces.PartialMedium
        "Medium model"
        annotation (choicesAllMatching=true,
          Documentation(info="<html>
        <u>Medium Model for the sensor. Make sure it is the same as for all lines the sensors input is connected.</u>
        </html>"));

      parameter Integer digits(min=0) = 1 "Number of displayed digits";
      parameter Boolean outputValue = false "Enable sensor-value output"
        annotation(Dialog(group="Output Value"));
      parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
        annotation(Dialog(group="Output Value", enable=outputValue));
      parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
        annotation(Dialog(tab="Initialization", enable=filter_output));
      parameter Real[Medium.nX] value_0(each unit="kg/kg") = Medium.X_default "Initial output state of sensor"
        annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
      parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));
      parameter Integer row(min=1, max=Medium.nX) = 1 "Row of mass fraction vector to display";

      Modelica.Blocks.Interfaces.RealOutput value_out[Medium.nX](each unit="kg/kg") = value if outputValue "Measured value [variable]"
        annotation (Placement(transformation(extent={{80,-20},{120,20}})));

      output Real value[Medium.nX](each unit="kg/kg") "Computed value of the selected Quantity";
      output Real display_value(unit="kg/kg") = value[row] "Row of the value vector to display";

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      Real direct_value[Medium.nX](each unit="kg/kg");

      function mfk = Chemical.Utilities.Functions.massFractionK (
                                                       redeclare package stateOfMatter =
              stateOfMatter);

    initial equation
      if filter_output and init==InitMode.steadyState then
        value= direct_value;
      elseif filter_output then
        value = value_0;
      end if;

    equation
      //OM fix
      //direct_value[1:Medium.nXi] = Medium.massFraction(state);

      for i in 1:Medium.nXi loop
        direct_value[i] = mfk(state, i);
      end for;

      if Medium.reducedX then
        direct_value[end] = 1-sum(direct_value[1:Medium.nXi]);
      end if;

      if filter_output then
        der(value) * TC = direct_value-value;
      else
        value = direct_value;
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-54,24},{66,-36}},
              lineColor={0,0,0},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{0,0},{0,-80}},
              color={158,66,200},
              thickness=0.5),
            Rectangle(
              extent={{-60,30},{60,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-60,30},{60,-30}},
              textColor={158,66,200},
              textString=DynamicSelect("value", String(display_value, format="1."+String(digits)+"f"))),
            Text(
              extent={{-26,22},{60,69}},
              textColor={175,175,175},
              textString="%row. mass-fraction"),
            Ellipse(
              extent={{-5,-75},{5,-85}},
              lineColor={158,66,200},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Sensor for measuring mass fraction X. Which row from X to display can be selected by the row parameter.</u>
<u>This sensor can be connected to a fluid stream without a junction.</u>
</html>"));
    end SingleSensorX;

    model SingleFlowSensor "Sensor for a selectable quantity associated with the massflow"
       extends Internal.PartialSensor;
      import Quantities=Chemical.Sensors.Internal.Types.MassFlowQuantities;
      import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

      parameter Quantities quantity "Quantity the sensor measures";
      parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimum Density" annotation (Dialog(tab="Advanced", group="Regularization"));
      parameter Boolean outputValue = false "Enable sensor-value output"
        annotation(Dialog(group="Output Value"));
      parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
        annotation(Dialog(group="Output Value", enable=outputValue));
      parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
        annotation(Dialog(tab="Initialization", enable=filter_output));
      parameter Real value_0(unit=Chemical.Sensors.Internal.getFlowUnit(quantity)) = 0 "Initial output state of sensor"
        annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
      parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));

      Modelica.Blocks.Interfaces.RealOutput value_out(unit=Chemical.Sensors.Internal.getFlowUnit(quantity)) = value if outputValue "Measured quantity [variable]"
        annotation (Placement(transformation(extent={{80,-20},{120,20}}),
            iconTransformation(extent={{80,-20},{120,20}})));

      output Real value(unit=Chemical.Sensors.Internal.getFlowUnit(quantity));

    protected
      Real direct_value(unit=Chemical.Sensors.Internal.getFlowUnit(quantity));

      function getQuantity = Chemical.Sensors.Internal.getFlowQuantity (
        redeclare package stateOfMatter = stateOfMatter) "Quantity compute function"
        annotation (Documentation(info="<html>
        <u>This function computes the selected quantity from state and massflow. rho_min is neddet for the computation of v. </u>
        </html>"));

    initial equation
      if filter_output and init==InitMode.steadyState then
        value= direct_value;
      elseif filter_output then
        value = value_0;
      end if;

    equation
      direct_value = getQuantity(state, rear.n_flow, quantity, rho_min);

      if filter_output then
        der(value) * TC = direct_value-value;
      else
        value = direct_value;
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-54,24},{66,-36}},
              lineColor={0,0,0},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{0,-26},{0,-80}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,-74},{6,-86}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Rectangle(
              extent={{-60,30},{60,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-60,30},{60,-30}},
              textColor={158,66,200},
              textString=DynamicSelect("value", String(value, format="1."+String(digits)+"f"))),
            Text(
              extent={{0,25},{60,75}},
              textColor={175,175,175},
              textString="%quantity")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>A undirected sensor measuring a selectable flow quantity associated with the massflow. For some quatities several units are available.</u>
</html>"));
    end SingleFlowSensor;

    model MultiSensor_Tpm "Undirected Sensor for Temperature, potential and mass-flow"
      extends Internal.PartialSensor;

      import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

      parameter Chemical.Sensors.Internal.Types.TemperatureUnit temperatureUnit="K" "Unit for the temperature output"
        annotation (
        Dialog(group="Units"),
        choicesAllMatching=true,
        Evaluate=true);
      parameter Chemical.Sensors.Internal.Types.ChemicalPotentialUnit potentialUnit="J/mol" "Unit for the potential output"
        annotation (
        Dialog(group="Units"),
        choicesAllMatching=true,
        Evaluate=true);
      parameter Chemical.Sensors.Internal.Types.MassFlowUnit massFlowUnit = "(kg/s)" "Unit for potential measurement and output"
        annotation(choicesAllMatching = true, Evaluate = true);
      parameter Boolean outputTemperature = false "Enable temperature output"
        annotation(Dialog(group="Output Value"));
      parameter Boolean outputChemicalPotential = false "Enable potential output"
        annotation(Dialog(group="Output Value"));
      parameter Boolean outputMassFlowRate = false "Enable massFlow output"
        annotation(Dialog(group="Output Value"));
      parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
        annotation(Dialog(group="Output Value", enable=(outputTemperature or outputChemicalPotential or outputMassFlowRate)));
      parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
        annotation(Dialog(tab="Initialization", enable=filter_output));
      parameter Real u_0(final quantity="ChemicalPotential", final unit=potentialUnit) = 0 "Initial output potential of sensor"
        annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
      parameter Real T_0(final quantity="ThermodynamicTemperature", final unit=temperatureUnit) = 0 "Initial output Temperature of sensor"
        annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
      parameter Real n_flow_0(final quantity="MassFlowRate", final unit=massFlowUnit) = 0 "Initial output massflow of sensor"
        annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
      parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant"
        annotation (Dialog(tab="Advanced", enable=(outputTemperature or outputChemicalPotential or outputMassFlowRate) and filter_output));

      Modelica.Blocks.Interfaces.RealOutput T_out(final quantity="ThermodynamicTemperature", final unit=temperatureUnit) = T if outputTemperature "Measured temperature [variable]"
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,60}),
            iconTransformation(extent={{80,60},{120,100}})));
      Modelica.Blocks.Interfaces.RealOutput u_out(final quantity="ChemicalPotential", final unit=potentialUnit) = u if outputChemicalPotential "Measured potential [variable]"
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,20}),
            iconTransformation(extent={{80,0},{120,40}})));
      Modelica.Blocks.Interfaces.RealOutput n_flow_out(unit="kg/s") = n_flow if outputMassFlowRate
        "Measured mass-flow [kg/s]"
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,-60}),
            iconTransformation(extent={{80,-60},{120,-20}})));

      output Real u(final quantity="ChemicalPotential", final unit=potentialUnit);
      output Real T(final quantity="ThermodynamicTemperature", final unit=temperatureUnit);
      output Real n_flow(final quantity="MassFlowRate", final unit=massFlowUnit);

    protected
      Real direct_p; // unit intentionally not set to avoid Warning
      Real direct_T; // unit intentionally not set to avoid Warning
      Real direct_n_flow; // unit intentionally not set to avoid Warning

    initial equation
      if filter_output and init==InitMode.steadyState then
        u=direct_p;
        T=direct_T;
        n_flow=direct_n_flow;
      elseif filter_output then
        u=u_0;
        T=T_0;
        n_flow=n_flow_0;
      end if;

    equation
      if temperatureUnit == "K" then
        direct_T =Medium.temperature(state);
      elseif temperatureUnit == "degC" then
        direct_T =Modelica.Units.Conversions.to_degC(Medium.temperature(state));
      end if;

      if potentialUnit == "J/mol" then
        direct_p =state.u;
      elseif potentialUnit == "bar" then
        direct_p =Modelica.Units.Conversions.to_bar(state);
      end if;

      if massFlowUnit == "(kg/s)" then
          direct_n_flow = rear.n_flow;
      elseif massFlowUnit == "(g/s)" then
          direct_n_flow = rear.n_flow*1000;
      end if;

      if filter_output then
        der(u) * TC = direct_p-u;
        der(T) * TC = direct_T-T;
        der(n_flow) * TC = direct_n_flow-n_flow;
      else
        u = direct_p;
        T = direct_T;
        n_flow = direct_n_flow;
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-54,94},{66,-66}},
              lineColor={0,0,0},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{-60,100},{60,-60}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{0,-60},{0,-80}}, color={0,0,0}),
            Ellipse(
              extent={{-6,-74},{6,-86}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{60,50}},
              textColor={158,66,200},
              textString=DynamicSelect("T", String(
                      T,
                      format="1."+String(digits)+"f"))),
            Text(
              extent={{-60,50},{60,0}},
              textColor={158,66,200},
              textString=DynamicSelect("u", String(
                      u,
                      format="1."+String(digits)+"f"))),
            Text(
              extent={{-60,0},{60,-50}},
              textColor={158,66,200},
              textString=DynamicSelect("m", String(
                      n_flow,
                      format="1."+String(digits)+"f"))),
            Text(
              extent={{-120,100},{-60,48}},
              textColor={175,175,175},
              textString="%temperatureUnit"),
            Text(
              extent={{-120,52},{-60,0}},
              textColor={175,175,175},
              textString="%potentialUnit"),
            Text(
              extent={{-120,0},{-60,-52}},
              textColor={175,175,175},
              textString="%massFlowUnit")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Undirected&nbsp;sensor&nbsp;for&nbsp;temperature,&nbsp;potential&nbsp;and&nbsp;mass-flow. Units can be selected.</u>
</html>"));
    end MultiSensor_Tpm;

    model UnidirectionalSensorAdapter "Adapter to connect a unidirectional sensor"
       extends Internal.PartialSensor;

      Chemical.Interfaces.Outlet outlet(redeclare package Medium=Medium)
        annotation (Placement(
            transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={0,40}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={0,-40})));

    equation
      outlet.r = rear.r;
      outlet.state = state;

      annotation (Icon(
        graphics={
          Line(
            points={{0,-80},{0,-40}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-5,-85},{5,-75}},
            lineColor={158,66,200},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)},
        coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,-20}})),
                                                      Diagram(
         coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,-20}})),
        Documentation(info="<html>
<u>A adapter to outputs the relevant state of the undirected flow, with r=0 at the outlet. It can be used to connect a unidirectional sensor to a undirected network.</u>
</html>"));
    end UnidirectionalSensorAdapter;

    package Tests "Tests package for the undirected sensor package"
    extends Modelica.Icons.ExamplesPackage;

      model TestSensors "Test for the undirected sensors"
        extends Modelica.Icons.Example;

        replaceable package Medium = Chemical.Media.myMedia.Water.StandardWater constrainedby
          Chemical.Media.myMedia.Interfaces.PartialTwoPhaseMedium
          "Medium model"
          annotation (Documentation(info="<html>
<u>Replaceable package with the medium model. Due to the vaporQuality sensor it must be a TwoPhaseMedium.</u>
</html>"));

        replaceable package Medium2 =
            Chemical.Media.myMedia.IdealGases.MixtureGases.CombustionAir                           constrainedby
          Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium model"
          annotation (Documentation(info="<html>
<u>Replaceable package with the medium model. Due to the vaporQuality sensor it must be a TwoPhaseMedium.</u>
</html>"));

        Chemical.Undirected.Processes.FlowResistance flowResistance(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.01,
          l=100,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-40,-8},{-20,12}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear(
          redeclare package stateOfMatter = stateOfMatter,
          T0_par=373.15,
          u0_par=100000) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-120,2})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par=303.15,
          u0_par=110000) annotation (Placement(transformation(extent={{84,-8},{104,12}})));
        inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=0.01) annotation (Placement(transformation(extent={{-130,22},{-110,42}})));
        Modelica.Blocks.Sources.Step step(
          height=-80000,
          offset=140000,
          startTime=5)
          annotation (Placement(transformation(extent={{120,-4},{108,8}})));
        MultiSensor_Tpm multiSensor_Tpm(redeclare package stateOfMatter =
              stateOfMatter,
          temperatureUnit="degC",
          potentialUnit="bar",
          outputTemperature=false)
          annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
        SingleSensorSelect singleSensorSelect(redeclare package stateOfMatter =
              stateOfMatter,
            quantity=Chemical.Sensors.Internal.Types.Quantities.u_bar)
          annotation (Placement(transformation(extent={{-10,0},{10,20}})));
        UnidirectionalSensorAdapter unidirectionalSensorAdapter(
          redeclare package stateOfMatter = stateOfMatter)
          annotation (Placement(transformation(extent={{20,0},{40,8}})));
        Chemical.Sensors.TwoPhaseSensorSelect sensor_vaporQuality1(
          redeclare package stateOfMatter = stateOfMatter, quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg)
          annotation (Placement(transformation(extent={{50,12},{70,32}})));
        SingleFlowSensor singleFlowSensor(redeclare package stateOfMatter =
              stateOfMatter,                                                               quantity=Chemical.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps)
          annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
        TwoPhaseSensorSelect sensor_vaporQuality2(
          redeclare package Medium2Phase = Medium,
          quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg,
          redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(extent={{50,0},{70,20}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance1(
          redeclare package stateOfMatter = stateOfMatter,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.01,
          l=100,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-40,52},{-20,72}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear1(
          redeclare package stateOfMatter = stateOfMatter,
          T0_par=373.15,
          u0_par=100000) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-120,62})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore1(
          redeclare package stateOfMatter = stateOfMatter,
          potentialFromInput=true,
          T0_par=303.15,
          u0_par=110000) annotation (Placement(transformation(extent={{84,52},{104,72}})));
        Modelica.Blocks.Sources.Step step1(
          height=-80000,
          offset=140000,
          startTime=5)
          annotation (Placement(transformation(extent={{120,56},{108,68}})));
        MultiSensor_Tpm multiSensor_Tpm1(
          redeclare package stateOfMatter = stateOfMatter,
          temperatureUnit="degC",
          potentialUnit="bar",
          outputTemperature=true,
          filter_output=true)
          annotation (Placement(transformation(extent={{-100,60},{-80,80}})));
        SingleSensorSelect singleSensorSelect1(
          redeclare package stateOfMatter = stateOfMatter,
          outputValue=true,
          quantity=Chemical.Sensors.Internal.Types.Quantities.u_bar,
          filter_output=true,
          init=Chemical.Sensors.Internal.Types.InitializationModelSensor.state,
          value_0=1) annotation (Placement(transformation(extent={{-10,60},{10,80}})));
        SingleFlowSensor singleFlowSensor1(
          redeclare package stateOfMatter = stateOfMatter,
          outputValue=true,
          quantity=Chemical.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps,
          filter_output=true) annotation (Placement(transformation(extent={{-70,60},{-50,80}})));
        TwoPhaseSensorSelect sensor_vaporQuality4(
          redeclare package Medium2Phase = Medium,
          outputValue=true,
          quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg,
          filter_output=true,
          redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(extent={{50,60},{70,80}})));
        UnidirectionalSensorAdapter unidirectionalSensorAdapter1(
          redeclare package stateOfMatter = stateOfMatter)
          annotation (Placement(transformation(extent={{20,64},{40,56}})));
        Chemical.Sensors.DifferenceTwoPhaseSensorSensorSelect differenceTwoPhaseSensorSensorSelect(
          redeclare package MediumA = Medium,
          redeclare package MediumB = Medium,
          quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.u_sat_Pa,
          outputValue=true,
          filter_output=true) annotation (Placement(transformation(extent={{50,52},{70,32}})));
        Chemical.Undirected.Boundaries.BoundaryRear boundary_rear2(
          redeclare package Medium = Medium2,
          T0_par=373.15,
          u0_par=200000,
          Xi0_par={0.2,0.8}) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-120,-30})));
        Chemical.Undirected.Boundaries.BoundaryFore boundary_fore2(
          redeclare package Medium = Medium2,
          potentialFromInput=false,
          T0_par=303.15,
          u0_par=100000) annotation (Placement(transformation(extent={{86,-40},{106,-20}})));
        Chemical.Undirected.Processes.FlowResistance flowResistance2(
          redeclare package Medium = Medium2,
          initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
          r=0.01,
          l=100,
          redeclare function pLoss =
              Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
              (                                                                                                        material=Chemical.Processes.Internal.Material.steel))
          annotation (Placement(transformation(extent={{-38,-40},{-18,-20}})));
        SingleSensorX singleSensorX(redeclare package Medium = Medium2) annotation (Placement(transformation(extent={{-100,-32},{-80,-12}})));
        SingleSensorX singleSensorX1(
          redeclare package Medium = Medium2,
          digits=2,
          outputValue=true,
          filter_output=true,
          row=2) annotation (Placement(transformation(extent={{-70,-32},{-50,-12}})));
        SensorState sensorState(redeclare package Medium = Medium2) annotation (Placement(transformation(extent={{20,-32},{40,-12}})));
      equation
        connect(step.y, boundary_fore.u0_var)
          annotation (Line(points={{107.4,2},{102,2},{102,8},{96,8}},
                                                         color={0,0,127}));
        connect(singleSensorSelect.rear, flowResistance.fore) annotation (Line(
            points={{-10,2},{-20,2}},
            color={158,66,200},
            thickness=0.5));
        connect(singleSensorSelect.fore, unidirectionalSensorAdapter.rear)
          annotation (Line(
            points={{10,2},{20,2}},
            color={158,66,200},
            thickness=0.5));
        connect(sensor_vaporQuality1.inlet, unidirectionalSensorAdapter.outlet)
          annotation (Line(
            points={{50,22},{30,22},{30,6}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm.fore, singleFlowSensor.rear) annotation (Line(
            points={{-80,2},{-70,2}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance.rear, singleFlowSensor.fore) annotation (Line(
            points={{-40,2},{-50,2}},
            color={158,66,200},
            thickness=0.5));
        connect(unidirectionalSensorAdapter.fore, sensor_vaporQuality2.rear)
          annotation (Line(
            points={{40,2},{50,2}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore.rear, sensor_vaporQuality2.fore)
          annotation (Line(
            points={{84,2},{70,2}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear.fore, multiSensor_Tpm.rear) annotation (Line(
            points={{-110,2},{-100,2}},
            color={158,66,200},
            thickness=0.5));
        connect(step1.y, boundary_fore1.u0_var) annotation (Line(points={{107.4,62},{102,62},{102,68},{96,68}},
                                                                                                color={0,0,127}));
        connect(singleSensorSelect1.rear, flowResistance1.fore)
          annotation (Line(
            points={{-10,62},{-20,62}},
            color={158,66,200},
            thickness=0.5));
        connect(multiSensor_Tpm1.fore, singleFlowSensor1.rear)
          annotation (Line(
            points={{-80,62},{-70,62}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance1.rear, singleFlowSensor1.fore)
          annotation (Line(
            points={{-40,62},{-50,62}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_fore1.rear, sensor_vaporQuality4.fore)
          annotation (Line(
            points={{84,62},{70,62}},
            color={158,66,200},
            thickness=0.5));
        connect(boundary_rear1.fore, multiSensor_Tpm1.rear)
          annotation (Line(
            points={{-110,62},{-100,62}},
            color={158,66,200},
            thickness=0.5));
        connect(singleSensorSelect1.fore, unidirectionalSensorAdapter1.rear)
          annotation (Line(
            points={{10,62},{20,62}},
            color={158,66,200},
            thickness=0.5));
        connect(unidirectionalSensorAdapter1.fore, sensor_vaporQuality4.rear)
          annotation (Line(
            points={{40,62},{50,62}},
            color={158,66,200},
            thickness=0.5));
        connect(differenceTwoPhaseSensorSensorSelect.inletA, unidirectionalSensorAdapter.outlet)
          annotation (Line(
            points={{50.4,38},{30,38},{30,6}},
            color={158,66,200},
            thickness=0.5));
        connect(differenceTwoPhaseSensorSensorSelect.inletB, unidirectionalSensorAdapter1.outlet)
          annotation (Line(
            points={{50.4,46},{30,46},{30,58}},
            color={158,66,200},
            thickness=0.5));
        connect(singleSensorX.rear, boundary_rear2.fore) annotation (Line(
            points={{-100,-30},{-110,-30}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance2.rear, singleSensorX1.fore) annotation (Line(
            points={{-38,-30},{-50,-30}},
            color={158,66,200},
            thickness=0.5));
        connect(singleSensorX1.rear, singleSensorX.fore) annotation (Line(
            points={{-70,-30},{-80,-30}},
            color={158,66,200},
            thickness=0.5));
        connect(flowResistance2.fore, sensorState.rear) annotation (Line(
            points={{-18,-30},{20,-30}},
            color={158,66,200},
            thickness=0.5));
        connect(sensorState.fore, boundary_fore2.rear) annotation (Line(
            points={{40,-30},{86,-30}},
            color={158,66,200},
            thickness=0.5));
        annotation (
          Icon(graphics,
               coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-60},{140,80}})),
          experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
          Documentation(info="<html>
<u>Test&nbsp;for&nbsp;the&nbsp;undirected&nbsp;sensors.</u>
<u>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
      end TestSensors;
      annotation (Documentation(info="<html>
<u>Tests package for the undirected sensor package.</u>
</html>"));
    end Tests;

    package Internal "Partials and functions"
      extends Modelica.Icons.InternalPackage;

      partial model PartialSensor "Partial undirected sensor"
        replaceable package Medium =
            Chemical.Media.myMedia.Interfaces.PartialMedium
          "Medium model" annotation (choicesAllMatching=true, Documentation(
              info="<html>
<u>Replaceable medium package for the sensor.</u>
</html>"));

        parameter Integer digits(min=0) = 1 "Number of displayed digits"
          annotation(Dialog(group="Sensor display"));
        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
          annotation (Dialog(tab="Advanced", group="Regularization"));

        Chemical.Undirected.Interfaces.Rear rear(redeclare package stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,-80})));
        Chemical.Undirected.Interfaces.Fore fore(redeclare package stateOfMatter =
              stateOfMatter)
          annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,-80})));

      /*  function regStepSt = Undirected.Internal.regStepState (
    redeclare package stateOfMatter = stateOfMatter) "RegStep function for a state"
    annotation (Documentation(info="<html>
<u><span style=\"font-family: Courier New;\">RegStep function for a state. The medium of the sensor is used and given to the function.</span></u>
</html>"));
*/
        Medium.ThermodynamicState state = Chemical.Interfaces.SubstanceState(u=u_reg,h= h_reg); //= regStepSt(rear.n_flow, rear.state_forwards, rear.state_rearwards, n_flow_reg);

      protected
        outer Chemical.DropOfCommons dropOfCommons;
       Modelica.Units.SI.ChemicalPotential u_reg=Chemical.Undirected.Internal.regStep(
                  rear.n_flow,
                  rear.state_forwards.u,
                  rear.state_rearwards.u,
                  n_flow_reg);
        Modelica.Units.SI.MolarEnthalpy h_reg=Chemical.Undirected.Internal.regStep(
                  rear.n_flow,
                  rear.state_forwards.h,
                  rear.state_rearwards.h,
                  n_flow_reg);





      equation
        for i in 1:Medium.nXi loop
          Xi_reg[i] = Chemical.Undirected.Internal.regStep(
                  rear.n_flow,
                  Xi_forwards[i],
                  Xi_rearwards[i],
                  n_flow_reg);
        end for;

        fore.state_forwards = rear.state_forwards;
        rear.state_rearwards = fore.state_rearwards;
        fore.r = rear.r;
        fore.n_flow + rear.n_flow = 0;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Line(
                points={{-100,-80},{100,-80}},
                color={158,66,200},
                thickness=0.5)}),
             Diagram(coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<u>Parent class of all undirected sensors.</u>
</html>"));
      end PartialSensor;
      annotation (Documentation(info="<html>
<u>Partials and functions needet for undirected sensors.</u>
</html>"));
    end Internal;
    annotation (Documentation(revisions="<html>
<u><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</u>

</html>",   info="<html>
<u>
Sensors package for undirected chemical simulation. For undirected flow, the massflow
must always be measured, since it determines if the forward or the backward flow is valid.
Therefore the fluid must flow through the sensor.
</u>
</html>"));
  end Sensors;

  package Internal "Internal helper functions and models for the undirected themofluid simulation."
    extends Modelica.Icons.InternalPackage;

    function regStep
      "Approximation of a general step, such that the characteristic is continuous and differentiable"
      extends Modelica.Icons.Function;
      input Real x "Abscissa value";
      input Real y1 "Ordinate value for x > 0";
      input Real y2 "Ordinate value for x < 0";
      input Real x_small(min=0) = 1e-5
        "Approximation of step for -x_small <= x <= x_small; x_small >= 0 required";
      output Real y "Ordinate value to approximate y = if x > 0 then y1 else y2";
    algorithm
      y := smooth(1, if x >  x_small then y1 else
                     if x < -x_small then y2 else
                     if x_small > 0 then (x/x_small)*((x/x_small)^2 - 3)*(y2-y1)/4 + (y1+y2)/2 else (y1+y2)/2);
      annotation(Documentation(revisions="<html>
<ul>
<li><em>April 29, 2008</em>
    by <a href=\"mailto:Martin.Otter@DLR.de\">Martin Otter</a>:<br>
    Designed and implemented.</li>
<li><em>August 12, 2008</em>
    by <a href=\"mailto:Michael.Sielemann@dlr.de\">Michael Sielemann</a>:<br>
    Minor modification to cover the limit case <code>x_small -> 0</code> without division by zero.</li>
</ul>
</html>",     info="<html>
<u>
This function is used to approximate the equation
</u>

<blockquote><pre>
y = <strong>if</strong> x &gt; 0 <strong>then</strong> y1 <strong>else</strong> y2;
</pre></blockquote>

<u>
by a smooth characteristic, so that the expression is continuous and differentiable:
</u>

<blockquote><pre>
y = <strong>smooth</strong>(1, <strong>if</strong> x &gt;  x_small <strong>then</strong> y1 <strong>else</strong>
              <strong>if</strong> x &lt; -x_small <strong>then</strong> y2 <strong>else</strong> f(y1, y2));
</pre></blockquote>

<u>
In the region -x_small &lt; x &lt; x_small a 2nd order polynomial is used
for a smooth transition from y1 to y2.
</u>
</html>"));
    end regStep;

    function regStepState "Apply regStep on State"
      extends Modelica.Icons.Function;

      input Modelica.Units.SI.MolarFlowRate n_flow;
      input Chemical.Interfaces.SubstanceState state_forwards;
      input Chemical.Interfaces.SubstanceState state_rearwards;
      input Modelica.Units.SI.MolarFlowRate n_flow_reg;

      output Chemical.Interfaces.SubstanceState state;

    protected
      Modelica.Units.SI.ChemicalPotential u;
      Modelica.Units.SI.MolarEnthalpy h;


    algorithm
      u := regStep(n_flow, state_forwards.u, state_rearwards.u, n_flow_reg);
      h := regStep(n_flow, state_forwards.h, state_rearwards.h, n_flow_reg);


      state := Chemical.Interfaces.SubstanceState(u=u, h=h);

      annotation (Documentation(info="<html>
<u>This function applies the regStep function to u,T and Xi of a state and creates and returns the resulting state.</u>
</html>"));
    end regStepState;
    annotation (Documentation(info="<html>
<u>Internal helper functions and models for the undirected themofluid simulation.</u>
</html>"));
  end Internal;
end Undirected;
