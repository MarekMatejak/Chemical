within Chemical;
package Interfaces "Chemical interfaces"
  import Chemical;

  connector Fore "Undirected connector outputting the forward state"

    Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
    flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

    Chemical.Interfaces.SubstanceStateOutput state_forwards "State of substance in forwards direction";
    Chemical.Interfaces.SubstanceStateInput state_rearwards "State of substance in rearwards direction";

    Chemical.Interfaces.SolutionStateOutput solution_forwards "State of solution to forward direction";
    Chemical.Interfaces.SolutionStateInput solution_rearwards "State of solution from rearward direction";

    Chemical.Interfaces.DefinitionOutput definition "Definition of substance";
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

    Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
    flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

    Chemical.Interfaces.SubstanceStateInput state_forwards "State of substance in forwards direction";
    Chemical.Interfaces.SubstanceStateOutput state_rearwards "State of substance in rearwards direction";

    Chemical.Interfaces.SolutionStateInput solution_forwards "State of solution from forward direction";
    Chemical.Interfaces.SolutionStateOutput solution_rearwards "State of solution to rearward direction";

    Chemical.Interfaces.DefinitionInput definition "Definition of substance";

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

  extends Modelica.Icons.InterfacesPackage;

  operator record Definition "Definition of a chemical substance or a chemical process"

    Chemical.Interfaces.DataRecord data "Data record of the substance or process";

     Boolean SelfClustering "default=false If true then the base molecules are binding together into clusters";
     Modelica.Units.SI.MolarEnthalpy SelfClustering_dH "Enthalpy of bond between two base molecules of substance at 25degC, 1 bar";
     Modelica.Units.SI.MolarEntropy SelfClustering_dS "Entropy of bond between two base molecules of substance at 25degC, 1 bar";


   encapsulated operator 'constructor'
      import Definition=Chemical.Interfaces.Definition;
      import ModelicaDataRecord=Modelica.Media.IdealGases.Common.DataRecord;
      import DataRecord=Chemical.Interfaces.DataRecord;
      import PhaseType=Chemical.Interfaces.Phase;
      //    constant Real R=1.380649e-23*6.02214076e23;
      //    constant Real T0=298.15 "Base temperature";
      //    constant Real p0=100000 "Base pressure";

      function fromDataRecord
        input DataRecord data "Mass based data record";
        input Boolean SelfClustering = false;
        input Real SelfClustering_dH = 0;
        input Real SelfClustering_dS = 0;
        output Definition result(
                  data=data,
                  SelfClustering=SelfClustering,
                  SelfClustering_dH=SelfClustering_dH,
                  SelfClustering_dS=SelfClustering_dS) "Molar based data record";
      algorithm
        annotation (Inline=true);
      end fromDataRecord;

      function fromFormationEnergies
        input Real MM=1 "Molar mass of the substance";
        input Real z=0 "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";
        input Real DfG=0 "Gibbs energy of formation of the substance at SATP conditions (298.15, 1bar)";
        input Real DfH=0 "Enthalpy of formation of the substance at SATP conditions (298.15, 1bar)";
        input Real S0=0 "Standart molar entropy of the substance at SATP conditions (298.15, 1bar)";
        input Real Cp=1 "Molar heat capacity of the substance at  SATP conditions (298.15, 1bar)";
        input PhaseType phase=PhaseType.Incompressible "State of matter";
        input Real Vm=if (phase == PhaseType.Gas) then ((1.380649e-23*6.02214076e23)*298.15)/100000 else 0.001*MM "Molar volume of the pure substance at SATP conditions (298.15, 1bar) (default fron non-gaseous is to reach density 1kg/L)";
        input Real gamma=1 "Activity coefficient of the substance";
        input Boolean SelfClustering = false;
        input Real SelfClustering_dH = 0;
        input Real SelfClustering_dS = 0;
        output Definition result(
                  data=DataRecord(
                    MM=MM,
                    Hf=DfH,
                    H0=DfH - 298.15*Cp,
                    alow={0,0,Cp/(1.380649e-23*6.02214076e23),0,0,0,0},
                    blow={(DfH-Cp*298.15)/(1.380649e-23*6.02214076e23),(S0-Cp*log(298.15))/(1.380649e-23*6.02214076e23),(((DfH-DfG)/298.15)-Cp*log(298.15))/(1.380649e-23*6.02214076e23)},
                    ahigh={0,0,Cp/(1.380649e-23*6.02214076e23),0,0,0,0},
                    bhigh={(DfH-Cp*298.15)/(1.380649e-23*6.02214076e23),(S0-Cp*log(298.15))/(1.380649e-23*6.02214076e23),(((DfH-DfG)/298.15)-Cp*log(298.15))/(1.380649e-23*6.02214076e23)},
                    z=z,
                    phase=phase,
                    VmBase=Vm/(1+log(gamma)),
                    VmExcess=Vm*log(gamma)/(1+log(gamma))),
                  SelfClustering=SelfClustering,
                  SelfClustering_dH=SelfClustering_dH,
                  SelfClustering_dS=SelfClustering_dS);
      algorithm
        annotation (Inline=true);
      end fromFormationEnergies;
   end 'constructor';


    encapsulated operator function '+'
      import Definition=Chemical.Interfaces.Definition;
      import DataRecord=Chemical.Interfaces.DataRecord;
      input Definition d1;
      input Definition d2;
      output Definition result " = d1 + d2";

    algorithm
      result :=Definition(
          data=DataRecord(
            MM=d1.data.MM + d2.data.MM,
            Hf=d1.data.Hf + d2.data.Hf,
            H0=d1.data.H0 + d2.data.H0,
            alow=d1.data.alow .+ d2.data.alow,
            blow=d1.data.blow .+ d2.data.blow,
            ahigh=d1.data.ahigh .+ d2.data.ahigh,
            bhigh=d1.data.bhigh .+ d2.data.bhigh,
            z=d1.data.z + d2.data.z,
            phase=d1.data.phase,
            VmBase=d1.data.VmBase + d2.data.VmBase,
            VmExcess=d1.data.VmExcess + d2.data.VmExcess));
          annotation (Inline=true);
    end '+';

    encapsulated operator '-'
      import Definition=Chemical.Interfaces.Definition;
      import DataRecord=Chemical.Interfaces.DataRecord;
     function negate
       input Definition d;
       output Definition result " = - d";
     algorithm
       result :=Definition(
          data=DataRecord(
            MM=-d.data.MM,
            Hf=-d.data.Hf,
            H0=-d.data.H0,
            alow=(-1) .* d.data.alow,
            blow=(-1) .* d.data.blow,
            ahigh=(-1) .* d.data.ahigh,
            bhigh=(-1) .* d.data.bhigh,
            z= -d.data.z,
            phase = d.data.phase,
            VmBase= -d.data.VmBase,
            VmExcess=d.data.VmExcess));

          annotation (Inline=true);
     end negate;

     function substract
      input Definition d1;
      input Definition d2;
      output Definition result " = d1 - d2";
     algorithm
      result :=Definition(
          data=DataRecord(
            MM=d1.data.MM - d2.data.MM,
            Hf=d1.data.Hf - d2.data.Hf,
            H0=d1.data.H0 - d2.data.H0,
            alow=d1.data.alow   .-  d2.data.alow,
            blow=d1.data.blow   .-  d2.data.blow,
            ahigh=d1.data.ahigh .-  d2.data.ahigh,
            bhigh=d1.data.bhigh .-  d2.data.bhigh,
            z=d1.data.z - d2.data.z,
            phase=d1.data.phase,
            VmBase=d1.data.VmBase - d2.data.VmBase,
            VmExcess=d1.data.VmExcess - d2.data.VmExcess));

          annotation (Inline=true);
     end substract;
    end '-';

    encapsulated operator '*'
      import Definition=Chemical.Interfaces.Definition;
      import DataRecord=Chemical.Interfaces.DataRecord;

    function scalar
      input Real n=1 "Stoichiometric coefficient";
      input Definition d;
      output Definition result " = n * d";
    algorithm
      result :=Definition(
          data=DataRecord(
            MM=n * d.data.MM,
            Hf=n * d.data.Hf,
            H0=n * d.data.H0,
            alow=n * d.data.alow,
            blow=n * d.data.blow,
            ahigh=n * d.data.ahigh,
            bhigh=n * d.data.bhigh,
            z=n * d.data.z,
            phase=d.data.phase,
            VmBase=n * d.data.VmBase,
            VmExcess=n * d.data.VmExcess));

          annotation (Inline=true);
    end scalar;

    function vector
      input Real[:] n "Stoichiometric coefficients";
      input Definition[:] d;
      output Definition result " = n * d";
    algorithm
      result :=Definition(
          data=DataRecord(
            MM=n * d.data.MM,
            Hf=n * d.data.Hf,
            H0=n * d.data.H0,
            alow= {sum({(n[i] * d[i].data.alow[j])  for i in 1:size(n,1)}) for j in 1:7},
            blow= {sum({(n[i] * d[i].data.blow[j])  for i in 1:size(n,1)}) for j in 1:3},
            ahigh={sum({(n[i] * d[i].data.ahigh[j]) for i in 1:size(n,1)}) for j in 1:7},
            bhigh={sum({(n[i] * d[i].data.bhigh[j]) for i in 1:size(n,1)}) for j in 1:3},
            z=n*d.data.z,
            phase=d[1].data.phase,
            VmBase=n*d.data.VmBase,
            VmExcess=n*d.data.VmExcess));
          annotation (Inline=true);
    end vector;
    end '*';

  end Definition;

 replaceable record SubstanceState "Set that defines a state of substance"
  extends Modelica.Icons.Record;

   Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential of the substance";
   Modelica.Units.SI.MolarEnthalpy h "Molar enthalpy of the substance";
 end SubstanceState;

 operator record SolutionState "Set that defines a state of solution"

   Modelica.Units.SI.Temperature T "Temperature of the solution";
   Modelica.Units.SI.Pressure p "Pressure of the solution";
   Modelica.Units.SI.ElectricPotential v "Electric potential in the solution";
   Modelica.Units.SI.AmountOfSubstance n "Amount of the solution";
   Modelica.Units.SI.Mass m "Mass of the solution";
   Modelica.Units.SI.Volume V "Volume of the solution";
   Modelica.Units.SI.Energy G "Free Gibbs energy of the solution";
   Modelica.Units.SI.ElectricCharge Q "Electric charge of the solution";
   Modelica.Units.SI.MoleFraction I "Mole fraction based ionic strength of the solution";

   encapsulated operator 'constructor'
     import Chemical.Interfaces.SolutionState;
     import Chemical.Interfaces.SolutionPort;
     import Chemical.Interfaces.Phase;
      //   constant Real R=1.380649e-23*6.02214076e23;
      //   constant Real T0=298.15 "Base temperature";
      //   constant Real p0=100000 "Base pressure";

     function fromValues
       input Phase phase "Phase of the chemical solution";
       input Real T=298.15 "Temperature of the solution";
       input Real p=100000 "Pressure of the solution";
       input Real v=0 "Electric potential in the solution";
       input Real n=1 "Amount of the solution";
       input Real m=1 "Mass of the solution";
       input Real V=if (phase==Phase.Gas) then n*(1.380649e-23*6.02214076e23)*T/p else 0.001 "Volume of the solution";
       input Real G=0 "Free Gibbs energy of the solution";
       input Real Q=0 "Electric charge of the solution";
       input Real I=0 "Mole fraction based ionic strength of the solution";
       output SolutionState result(T=T,p=p,v=v,n=n,m=m,V=V,G=G,Q=Q,I=I);
     algorithm
       annotation(Inline = true);
     end fromValues;


   end 'constructor';

   encapsulated operator function '=='
     import Chemical.Interfaces.SolutionState;
     input SolutionState s1;
     input SolutionState s2;
     output Boolean result "= s1 == s2";
   algorithm
      result := s1.T == s2.T and s1.p == s2.p and s1.v == s2.v and s1.n == s2.n and s1.m == s2.m and s1.V == s2.V and s1.G == s2.G and s1.Q == s2.Q and s1.I == s2.I;
      annotation(Inline = true);
   end '==';
 end SolutionState;

 function processData "Process changes of Gibbs energy, enthalpy, volume and heat capacity (products - reactants)"
   import Chemical;

      extends Modelica.Icons.Function;
  input Real K "Process dissociation constant (mole-fraction based) at 25°C,1bar";
  input Modelica.Units.SI.MolarEnergy dH=0 "Process molar enthalpy change at 25°C,1bar";
  input Modelica.Units.SI.MolarHeatCapacity dCp=0 "Process molar heat capacity change at 25°C,1bar";
  input Modelica.Units.SI.MolarVolume dVm=0 "Process molar volume change at 25°C,1bar";
   input Chemical.Interfaces.Phase phase=Chemical.Interfaces.Phase.Incompressible "State of matter";
   output Chemical.Interfaces.Definition processData "Data record of process changes";
 algorithm
    processData := Chemical.Interfaces.Definition(
       MM=0,
       z=0,
       DfG=-Modelica.Constants.R*Chemical.Interfaces.Properties.T0*log(K),
       DfH=dH,
       Cp=dCp,
       Vm=dVm,
       phase=phase,
       gamma=1);

 end processData;

  connector SolutionPort
  "Only for connecting the one solution their substances. Please, do not use it in different way."

    //enthalpy
  Modelica.Units.SI.Temperature T "Temperature of the solution";
  flow Modelica.Units.SI.EnthalpyFlowRate dH
    "Internal enthalpy change of the solution";

    //pressure
  Modelica.Units.SI.Pressure p "Pressure of the solution";
  flow Modelica.Units.SI.VolumeFlowRate dV
    "Volume change of the solution";

    //electric port
  Modelica.Units.SI.ElectricPotential v
    "Electric potential in the solution";
  flow Modelica.Units.SI.ElectricCurrent i "Change of electric charge";

    //Extensive properties of the solution:

    // The extensive quantities here have not the real physical flows.
    // They hack the Kirchhof's flow equation to be counted as the sum from all connected substances in the solution.

    //amount of substances
  Modelica.Units.SI.AmountOfSubstance n "Amount of the solution";
  flow Modelica.Units.SI.AmountOfSubstance nj
    "Amount of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //mass of substances
  Modelica.Units.SI.Mass m "Mass of the solution";
  flow Modelica.Units.SI.Mass mj
    "Mass of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //volume of substances
  Modelica.Units.SI.Volume V "Volume of the solution";
  flow Modelica.Units.SI.Volume Vj
    "Volume of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //Gibbs energy of substances
  Modelica.Units.SI.Energy G "Free Gibbs energy of the solution";
  flow Modelica.Units.SI.Energy Gj
    "Gibbs energy of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //electric charge of the substance
  Modelica.Units.SI.ElectricCharge Q "Electric charge of the solution";
  flow Modelica.Units.SI.ElectricCharge Qj
    "Electric charge of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //ionic strength of substances
  Modelica.Units.SI.MoleFraction I
    "Mole fraction based ionic strength of the solution";
  flow Modelica.Units.SI.MoleFraction Ij
    "Mole-fraction based ionic strength of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

  /*  //suport for structural properties
  replaceable package stateOfMatter = StateOfMatter  constrainedby StateOfMatter
  "Substance model to translate data into substance properties"
     annotation (choicesAllMatching = true);*/

    annotation (
    defaultComponentName="solution",
    Documentation(revisions="<html>
<p><i>2015-2016</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Solution port integrates all substances of the solution:</p>
<p>Such as if there are connected together with electric port, thermal port and with port composed with the amont of substance and molar change of substance.</p>
</html>"),
         Icon(graphics={            Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={127,127,0},
            fillColor={127,127,0},
            fillPattern=FillPattern.Solid)}),
    Diagram(graphics={
     Text(extent={{-160,110},{40,50}},   lineColor={127,127,0},    textString = "%name",
          fillColor={127,127,0},
          fillPattern=FillPattern.Solid),
                    Rectangle(
            extent={{-40,40},{40,-40}},
            lineColor={127,127,0},
            fillColor={127,127,0},
            fillPattern=FillPattern.Solid,
            lineThickness=1)}));
  end SolutionPort;

  type Phase                 = enumeration(
    Gas
      "Gaseous phase",
    Liquid
      "Liquid phase",
    Solid
      "Liquid phase",
    Aqueous
      "Dissolved in water",
    Incompressible
      "Incompressible liquid or solid phase")
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(
        coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>
Phase of substance or solution. It is possible to include any phase - as default it behave as incompressible. 
To change its behavior it is necessary to modify Property functions.
</p>
</html>"));

  model ProcessProperties "Properties of the chemical process"
    import Chemical;

    Chemical.Interfaces.DefinitionInput definition "Definition of the process";

    Chemical.Interfaces.SolutionStateInput solutionState "State of the solution";

    Modelica.Units.SI.ChemicalPotential dG "Gibbs energy change during the process";

    Modelica.Units.SI.MolarEnthalpy dH "Molar enthalpy change during the process";

    Modelica.Units.SI.Temperature dlnK_per_dinvT(displayUnit="K") "Temperature dependence coefficient";

    Modelica.Units.SI.MolarEntropy dS "Molar entropy change during the process";

    Modelica.Units.SI.MolarHeatCapacity dCp "Molar heat capacity change during the process";

    Modelica.Units.SI.MolarVolume dVm "Molar volume change during the process";

    Modelica.Units.SI.MolarVolume dVmExcess "Molar volume excess change during the process";

    Real K "Dissociation constant (mole-fraction based)";

    outer Modelica.Fluid.System system "System wide properties";

  equation
   assert(abs(definition.data.MM) < 1e-5, "Process should not change the mass");
    assert(abs(Chemical.Interfaces.Properties.chargeNumberOfIon(definition, solutionState)) < Modelica.Constants.eps, "Process should not change the charge");

    dH = Chemical.Interfaces.Properties.molarEnthalpy(definition, solutionState);
    dlnK_per_dinvT = -dH/Chemical.Interfaces.Properties.R;

    dG = Chemical.Interfaces.Properties.chemicalPotentialPure(definition, solutionState);
   dG = dH - solutionState.T*dS;//  dS = molarEntropy(definition,dG,solutionState);

    K = exp(-dG/(Chemical.Interfaces.Properties.R*solutionState.T));

    dCp = Chemical.Interfaces.Properties.molarHeatCapacityCp(definition, solutionState);
    dVm = Chemical.Interfaces.Properties.molarVolume(definition, solutionState);
    dVmExcess = Chemical.Interfaces.Properties.molarVolumeExcess(definition, solutionState);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ProcessProperties;

  package Properties "Calculation of chemical properties from definition and solution state"

    constant Real R=1.380649e-23*6.02214076e23;
    constant Real T0=298.15 "Base temperature";
    constant Real p0=100000 "Base pressure";

    model SubstanceProperties "Properties of the substance"
      import Chemical;

      parameter Boolean FixedDefinition "definition==definitionParam";

      parameter Chemical.Interfaces.Definition definitionParam "used only if FixedDefinition to help initialization";

      parameter Modelica.Units.SI.Mass m_start "Start value for mass of the substance";

      parameter Boolean SolutionObserverOnly = false "True if substance does not affect the solution";

      Chemical.Interfaces.DefinitionInput definition "Definition of the substance";

      Chemical.Interfaces.SolutionStateInput solutionState "State of the solution";

      InputAmountOfSubstance amountOfBaseMolecules
        "Amount of base molecules inside all clusters in compartment";

      InputMolarFlowRate n_flow "Molar change of the substance";

      InputHeatFlowRate h_flow "Substance enthaply change";

      Modelica.Units.SI.MoleFraction x "Mole fraction of the base molecule of the substance";

      Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
            "Molar concentration of particles";

      Modelica.Units.SI.MassConcentration M(displayUnit="mg/l")
            "Mass concentration";

      Modelica.Units.SI.Molality b(displayUnit="mmol/kg")
            "Molality";

      Modelica.Units.SI.MassFraction X
            "Mass fraction";

      Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential of the substance";

      Modelica.Units.SI.ActivityOfSolute a
      "Activity of the substance (mole-fraction based)";

      Modelica.Units.SI.ActivityCoefficient gamma
      "Activity coefficient of the substance";

      Modelica.Units.SI.ChargeNumberOfIon z "Charge number of ion";

      Modelica.Units.SI.MolarEnthalpy h
      "Molar enthalpy of the substance";

      Modelica.Units.SI.MolarEnthalpy hPure
      "Molar enthalpy of the pure substance";

      Modelica.Units.SI.MolarEntropy sPure
      "Molar entropy of the pure substance";

      Modelica.Units.SI.ChemicalPotential u0
      "Chemical potential of the pure substance";

      Modelica.Units.SI.ChemicalPotential uPure
      "Electro-Chemical potential of the pure substance";

      Modelica.Units.SI.MolarVolume Vm
      "Molar volume of the substance";

      Modelica.Units.SI.MolarMass MM
      "Molar mass of the substance";

      Modelica.Units.SI.MolarVolume VmPure
      "Molar volume of the pure substance";

      Modelica.Units.SI.MolarVolume VmExcess
      "Molar volume excess of the substance in solution (typically it is negative as can be negative)";

      Modelica.Units.SI.SpecificVolume Vs "Specific volume";

      Modelica.Units.SI.Density rho "Density";

     // Local connector definition, used for equation balancing check

      connector InputMolarFlowRate = input Modelica.Units.SI.MolarFlowRate
        "Molar flow rate as input signal connector";
      connector InputHeatFlowRate = input Modelica.Units.SI.HeatFlowRate
        "Heat flow rate as input signal connector";

      connector InputSubstanceData = input Chemical.Interfaces.Definition
        "Substance definition as input signal connector";
      connector InputTemperature = input Modelica.Units.SI.Temperature
        "Temperature as input signal connector";
      connector InputPressure = input Modelica.Units.SI.Pressure
        "Pressure as input signal connector";
      connector InputElectricPotential = input
          Modelica.Units.SI.ElectricPotential
        "Electric potential as input signal connector";
      connector InputMoleFraction = input Modelica.Units.SI.MoleFraction
        "Mole fraction as input signal connector";
      connector InputAmountOfSubstance = input
          Modelica.Units.SI.AmountOfSubstance
        "Amount of substance as input signal connector";

      Modelica.Units.SI.AmountOfSubstance amountOfFreeMolecule(start=
           m_start*specificAmountOfFreeBaseMolecule(
                                       definitionParam,
                                       Chemical.Interfaces.SolutionState(
                                       phase=Chemical.Interfaces.Phase.Incompressible,
                                       T=system.T_ambient,
                                       p=system.p_ambient)))
        "Amount of free molecules not included inside any clusters in compartment";

      Modelica.Units.SI.AmountOfSubstance amountOfParticles(start=
           m_start*specificAmountOfParticles(
                                       definitionParam,
                                       Chemical.Interfaces.SolutionState(
                                       phase=Chemical.Interfaces.Phase.Incompressible,
                                       T=system.T_ambient,
                                       p=system.p_ambient)))
        "Amount of particles/clusters in compartment";

      Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
          Modelica.Constants.R*solutionState.T))
        "Dissociation constant of hydrogen bond between base molecules";

      Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
          selfClusteringBondEnthalpy(definition)
        - solutionState.T * selfClusteringBondEntropy(definition)
        "Gibbs energy of hydrogen bond between H2O molecules";

      Modelica.Units.SI.AmountOfSubstance amountOfBonds
        "Amount of hydrogen bonds between molecules in compartment";

      outer Modelica.Fluid.System system "System wide properties";

      Real i,dH,dV,nj,mj,Vj,Gj,Qj,Ij;

    equation
     // assert(x > 0, "Molar fraction must be positive");

     gamma = activityCoefficient(definition,solutionState);
     z = chargeNumberOfIon(definition,solutionState);

     MM = 1/specificAmountOfParticles(definition,solutionState);
     h = molarEnthalpy(definition,solutionState);
     sPure = molarEntropyPure(definition,solutionState);
     hPure = uPure + solutionState.T*sPure;
     u0 = chemicalPotentialPure(definition,solutionState);
     uPure = electroChemicalPotentialPure(definition,solutionState);
     Vm = molarVolume(definition,solutionState);
     VmPure =Chemical.Interfaces.Properties.molarVolumeBase(definition, solutionState);
     VmExcess = molarVolumeExcess(definition,solutionState);
     //  molarHeatCapacityCp = smolarHeatCapacityCp(definition,solutionState);

     //activity of the substance
     a = gamma*x;

     //electro-chemical potential of the substance in the solution
     u = chemicalPotentialPure(definition,solutionState)
       + (Modelica.Constants.R*solutionState.T)*log(a)
       + z*Modelica.Constants.F*solutionState.v;

      // during initialization it does not take value from parameter definition if not useInlet,
      // so instead of just "selfClustering(definition)" it must be written
      if (FixedDefinition and selfClustering(definitionParam)) or selfClustering(definition) then

        //Liquid cluster theory - equilibrium:
        //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules
        //amountOfParticles/solutionState.n = x/(1-K*x);                //sum(x[i])
        //amountOfBaseMolecules/solutionState.n = x/((1-K*x)^2);            //sum(i*x[i])
        //amountOfHydrogenBonds/solutionState.n = x*x*K/((1-K*x)^2);   //sum((i-1)*x[i])

        amountOfParticles*(1 - SelfClustering_K*x) = amountOfFreeMolecule;

        //Calculation of "abs(amountOfBaseMolecules*(1 - SelfClustering_K*x)) = amountOfParticles":
        x = ((2*SelfClustering_K+solutionState.n/amountOfBaseMolecules) - sqrt((4*SelfClustering_K*solutionState.n/amountOfBaseMolecules)+(solutionState.n/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        amountOfBonds = amountOfBaseMolecules*x*SelfClustering_K;

        //TODO: may be the volume of the same number of free water molecules is different as volume of the same number of water molecules in cluster ..
        //TODO: more precise calculation of other properties

       //der(enthalpy) = dH + n_flow*actualStream(port_a.h_outflow);
       //enthalpy = h*amountOfBaseMolecules + amountOfAdditionalBonds*bondEnthalpy;
        dH = der(h)*
          amountOfBaseMolecules + n_flow*h - h_flow +
          selfClusteringBondEnthalpy(definition)*der(amountOfBonds)
                        "heat transfer from other substances in solution [J/s]";

        Gj =amountOfBaseMolecules*u + amountOfBonds*SelfClustering_dG
                        "Gibbs energy of the substance";

      else

        amountOfParticles = amountOfFreeMolecule;
        amountOfBaseMolecules = amountOfFreeMolecule;
        amountOfBonds = 0;

        //der(enthalpy) = dH + n_flow*actualStream(port_a.h_outflow);
        //enthalpy = h*amountOfBaseMolecules;
        dH =  der(h)*amountOfBaseMolecules + n_flow*h - h_flow
                  "heat transfer from other substances in solution [J/s]";

        Gj = amountOfBaseMolecules*u "Gibbs energy of the substance [J]";

      end if;

      x = amountOfFreeMolecule/solutionState.n "mole fraction of base molecule [mol/mol]";

      c = amountOfParticles/solutionState.V "concentration [mol/m3]";

      M = (amountOfParticles/specificAmountOfParticles(definition,solutionState))/solutionState.V "mass concentration [kg/m3]";

      b = amountOfParticles/solutionState.m "molality [mol/kg]";

      X  = (amountOfParticles/specificAmountOfParticles(definition,solutionState))/solutionState.m "mass fraction [kg/kg]";

      Vs = specificVolume(definition,solutionState);

      rho = density(definition,solutionState);

      if SolutionObserverOnly then
        i = 0;
        dV = 0;
        nj = 0;
        mj = 0;
        Vj = 0;
        Qj = 0;
        Ij = 0;
      else
      //solution flows
        i = Modelica.Constants.F*z*n_flow +
          Modelica.Constants.F*der(z)*amountOfBaseMolecules "change of sunstance charge [A]";
        dV = Vm*n_flow + der(Vm)*amountOfBaseMolecules "change of substance volume [m3/s]";

      //extensive properties
        nj = amountOfParticles;
        mj = amountOfBaseMolecules*molarMassOfBaseMolecule(definition);
        Vj = amountOfBaseMolecules*Vm;
        Qj = Modelica.Constants.F*amountOfBaseMolecules*z;
        Ij = (1/2)*(amountOfBaseMolecules*z^2);
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SubstanceProperties;

   function firstProductDefinition
       "Return formation definition of the first product of chemical process"
     import Chemical;

        extends Modelica.Icons.Function;
    //input
    input Modelica.Units.SI.StoichiometricNumber s[:] "Stoichiometric reaction coefficient for substrates [nS]";
    input Modelica.Units.SI.StoichiometricNumber p[:] "Stoichiometric reaction coefficient for products [nP]";
     input Chemical.Interfaces.Definition process "Data record of process changes";
     input Chemical.Interfaces.Definition substrates[:] "Substrates definitions [nS]";
      input Chemical.Interfaces.Definition products[:] "Products definitions [nP]";
     output Chemical.Interfaces.Definition firstProductDefinition "Definition of the first product in process";
   protected
     Chemical.Interfaces.Definition pd[size(products,1)-1];
   algorithm
        for i in 1:size(products,1)-1 loop
          pd[i].data := products[i+1].data;
        end for;
        firstProductDefinition := (1/p[1]) * (s*substrates + process - p[2:end]*pd);
         /*Definition(
       MM = ((p[2:end]*products.data.MM) - (s*substrates.MolarWeight))/p[1],
       z = ((p[2:end]*products.z) - (s*substrates.z))/p[1],
       DfG = ((p[2:end]*products.DfG) - (s*substrates.DfG) - process.DfG)/p[1],
       DfH = ((p[2:end]*products.DfH) - (s*substrates.DfH) - process.DfH)/p[1],
       gamma = 1,
       Cp = ((p[2:end]*products.Cp) - (s*substrates.Cp) - process.Cp)/p[1],
       SelfClustering = 0,
       SelfClustering_dH = 0,
       SelfClustering_dS = 0,
       Vs = ((p[2:end]*products.Vs) - (s*substrates.Vs) - process.Vs)/p[1]);*/
   end firstProductDefinition;

   function moleFraction "Return mole fraction from potential, definition and solution"
     import Chemical;
      extends Modelica.Icons.Function;
      input Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential";
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      output Modelica.Media.Interfaces.Types.MoleFraction moleFraction "Mole fraction of the substance";
   algorithm
       moleFraction := exp((u - chemicalPotentialPure(definition,solution)-definition.data.z*Modelica.Constants.F*solution.v)/(Modelica.Constants.R*solution.T))/activityCoefficient(definition,solution);
   end moleFraction;

   function activity "Return activity from potential, definition and solution"
     import Chemical;
      extends Modelica.Icons.Function;
      input Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential";
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      output Real activity "Activity of the substance";
   algorithm
       activity := exp((u - chemicalPotentialPure(definition,solution)-definition.data.z*Modelica.Constants.F*solution.v)/(Modelica.Constants.R*solution.T));
   end activity;

   function activityCoefficient
    "Return activity coefficient of the substance in the solution"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      output Real activityCoefficient "Activity Coefficient";
   algorithm
       activityCoefficient := if definition.data.phase==Chemical.Interfaces.Phase.Gas then 1 else exp(definition.data.VmExcess/definition.data.VmBase);
   end activityCoefficient;

   function chargeNumberOfIon
    "Return charge number of the substance in the solution"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      output Modelica.Units.SI.ChargeNumberOfIon chargeNumberOfIon
      "Charge number of ion";
   algorithm
       chargeNumberOfIon := definition.data.z;
   end chargeNumberOfIon;

   function molarEnthalpyElectroneutral
    "Molar enthalpy of the substance in electroneutral solution"
     import Chemical;
      extends Modelica.Icons.Function;
      import Modelica.Math;
     input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";

      output Modelica.Units.SI.MolarEnthalpy molarEnthalpyElectroneutral
      "Molar enthalpy";
   algorithm
       //Molar enthalpy:
       // - temperature and pressure shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - dp*Vm)
       //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - (dp/dT)*Vm. As a result dH = dT*Cv - dp*Vm for incompressible substances.

       molarEnthalpyElectroneutral := smooth(0,(if solution.T < Tlimit then R*((-definition.data.alow[1] + solution.T*(definition.data.blow[1] +
       definition.data.alow[2]*Math.log(solution.T) + solution.T*(1.*definition.data.alow[3] + solution.T*(0.5*definition.data.alow[4] +
       solution.T*(1/3*definition.data.alow[5] + solution.T*(0.25*definition.data.alow[6] + 0.2*definition.data.alow[7]*solution.T))))))
       /solution.T) else R*((-definition.data.ahigh[1] + solution.T*(definition.data.bhigh[1] + definition.data.ahigh[2]*
       Math.log(solution.T) + solution.T*(1.*definition.data.ahigh[3] + solution.T*(0.5*definition.data.ahigh[4] + solution.T*(1/3*definition.data.ahigh[5] +
       solution.T*(0.25*definition.data.ahigh[6] + 0.2*definition.data.ahigh[7]*solution.T))))))/solution.T)));


   end molarEnthalpyElectroneutral;

   function molarEnthalpy
    "Molar enthalpy of the substance with electric potential dependence"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      output Modelica.Units.SI.MolarEnthalpy molarEnthalpy "Molar enthalpy";
   algorithm
      molarEnthalpy := molarEnthalpyElectroneutral(definition,solution) +
           Modelica.Constants.F*chargeNumberOfIon(definition,solution)*solution.v;
      annotation (Inline=true, smoothOrder=2);
   end molarEnthalpy;

   function molarEntropyPure
    "Molar entropy of the pure substance"
     import Chemical;
      extends Modelica.Icons.Function;
      import Modelica.Math;
     input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.MolarEntropy molarEntropyPure
      "Molar entropy of the pure substance";

   algorithm
       //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
       //+ substanceData.Cv*log(solution.T/298.15);

       //Molar entropy shift:
       // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = solution.T*dS (small amount of added heat energy)
       // - pressure shift: with constant molar volume at constant temperature Vm*dP = -solution.T*dS (small amount of work)
       molarEntropyPure :=(if solution.T < Tlimit
        then
        R*(definition.data.blow[3] - 0.5*definition.data.alow[1]/(solution.T*solution.T)
         - definition.data.alow[2]/solution.T + definition.data.alow[3]*Math.log(solution.T)
         + solution.T*(definition.data.alow[4] + solution.T*(0.5*definition.data.alow[5]
         + solution.T*(1/3*definition.data.alow[6] + 0.25*definition.data.alow[7]*solution.T))))
        else
        R*(definition.data.bhigh[3] - 0.5*definition.data.ahigh[1]/(solution.T*solution.T)
         - definition.data.ahigh[2]/solution.T + definition.data.ahigh[3]*Math.log(solution.T)
         + solution.T*(definition.data.ahigh[4] + solution.T*(0.5*definition.data.ahigh[5]
         + solution.T*(1/3*definition.data.ahigh[6] + 0.25*definition.data.ahigh[7]*solution.T)))))
        + (if definition.data.phase == Chemical.Interfaces.Phase.Gas
             then -R*log(solution.p/100000)
             else -(Chemical.Interfaces.Properties.molarVolumeBase(definition, solution)/solution.T)*(solution.p - 100000));


       //For example at triple point of water should be solution.T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
       //As definition.data: http://www1.lsbu.ac.uk/water/water_phase_diagram.html
       //At solution.T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K

   end molarEntropyPure;

    function molarEntropy "Molar entropy of the substance in the solution"
      import Chemical;
          extends Modelica.Icons.Function;
      input Modelica.Units.SI.ChemicalPotential u
      "Electro-chemical potential of the substance";
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";

      output Modelica.Units.SI.MolarEntropy molarEntropy "Molar entropy";
    algorithm
        molarEntropy :=  (u - molarEnthalpy(definition,solution))/solution.T;
    end molarEntropy;

   function chemicalPotentialPure "Chemical potential of the pure substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.ChemicalPotential chemicalPotentialPure
      "Base chemical potential";
   algorithm
       chemicalPotentialPure :=  molarEnthalpyElectroneutral(definition,solution) - solution.T*molarEntropyPure(definition,solution);
   end chemicalPotentialPure;

   function electroChemicalPotentialPure
    "Electro-chemical potential of the pure substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.ChemicalPotential
      electroChemicalPotentialPure "Base electro-chemical potential";
   algorithm
    electroChemicalPotentialPure := chemicalPotentialPure(
         definition,
         solution) + Modelica.Constants.F*chargeNumberOfIon(definition,solution)*solution.v;
   end electroChemicalPotentialPure;

   function molarVolumeBase "Base molar volume of the substance"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.MolarVolume molarVolumePure "Molar volume";
   algorithm
     molarVolumePure :=(if definition.data.phase == Chemical.Interfaces.Phase.Gas then R*solution.T/solution.p else definition.data.VmBase);

   end molarVolumeBase;

   function molarVolumeExcess
    "Excess molar volume of the substance in the solution"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.MolarVolume molarVolumeExcess
      "Excess molar volume of the substance in the solution";
   algorithm
     //zero if activityCoefficient==1
      molarVolumeExcess :=(if definition.data.phase == Chemical.Interfaces.Phase.Gas then 0 else definition.data.VmExcess);

      annotation (Inline=true, smoothOrder=2);
   end molarVolumeExcess;

   function molarVolume "Molar volume of the substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.MolarVolume molarVolume "Molar volume";
   algorithm
    molarVolume :=Chemical.Interfaces.Properties.molarVolumeBase(definition, solution) + molarVolumeExcess(definition, solution);
      annotation (Inline=true, smoothOrder=2);
   end molarVolume;

   function molarHeatCapacityCp
    "Molar heat capacity at constant pressure"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.MolarHeatCapacity molarHeatCapacityCp
      "Molar heat capacity at constant pressure";
   algorithm
     molarHeatCapacityCp := smooth(0,if solution.T < Tlimit then R*(1/(solution.T*solution.T)*(definition.data.alow[1] + solution.T*(
       definition.data.alow[2] + solution.T*(1.*definition.data.alow[3] + solution.T*(definition.data.alow[4] + solution.T*(definition.data.alow[5] + solution.T
       *(definition.data.alow[6] + definition.data.alow[7]*solution.T))))))) else R*(1/(solution.T*solution.T)*(definition.data.ahigh[1]
        + solution.T*(definition.data.ahigh[2] + solution.T*(1.*definition.data.ahigh[3] + solution.T*(definition.data.ahigh[4] + solution.T*(definition.data.
       ahigh[5] + solution.T*(definition.data.ahigh[6] + definition.data.ahigh[7]*solution.T))))))));
   end molarHeatCapacityCp;

   function molarMassOfBaseMolecule
      "Molar mass of base molecule of the substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Data record of substance";
    output Modelica.Units.SI.MolarMass molarMass "Molar mass";
   algorithm
     molarMass := definition.data.MM;
   end molarMassOfBaseMolecule;

   function selfClustering "returns true if substance molecules are joining together to clusters"
     import Chemical;
       extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Data record of substance";
          output Boolean selfClustering;
   algorithm
     selfClustering:=definition.SelfClustering;
   end selfClustering;

   function selfClusteringBondEnthalpy
    "Enthalpy of joining two base molecules of the substance together to cluster"
     import Chemical;
       extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Data record of substance";
    output Modelica.Units.SI.MolarEnthalpy selfClusteringEnthalpy;
   algorithm
     selfClusteringEnthalpy:=definition.SelfClustering_dH;
   end selfClusteringBondEnthalpy;

   function selfClusteringBondEntropy
    "Entropy of joining two base molecules of the substance together to cluster"
     import Chemical;
       extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Data record of substance";
    output Modelica.Units.SI.MolarEntropy selfClusteringEntropy;
   algorithm
     selfClusteringEntropy:=definition.SelfClustering_dS;
   end selfClusteringBondEntropy;

   function selfClusteringBondVolume
     import Chemical;
       extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Data record of substance";
    output Modelica.Units.SI.MolarVolume selfClusteringBondVolume;
   algorithm
     selfClusteringBondVolume:=0;
   end selfClusteringBondVolume;

   function selfClusteringBondHeatCapacityCp
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Data record of substance";
    output Modelica.Units.SI.MolarHeatCapacity selfClusteringBondHeatCapacityCp;
   algorithm
     selfClusteringBondHeatCapacityCp:=0;
   end selfClusteringBondHeatCapacityCp;

    function specificAmountOfParticles
      "Amount of particles per mass of the substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfSubstance(unit="mol/kg")
        "Amount of substance particles per its mass";

    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
      Real amountOfBaseMolecules,amountOfFreeMolecule,amountOfParticles;
      Real x;
    algorithm
      if not selfClustering(definition) then
        specificAmountOfSubstance := 1/molarMassOfBaseMolecule(definition);
      else
        SelfClustering_dG :=selfClusteringBondEnthalpy(definition) - solution.T*
          selfClusteringBondEntropy(definition);

        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));

        amountOfBaseMolecules:=mass/molarMassOfBaseMolecule(definition);
        x:=((2*SelfClustering_K+nSolution/amountOfBaseMolecules) -
         sqrt((4*SelfClustering_K*nSolution/amountOfBaseMolecules)+
         (nSolution/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        amountOfFreeMolecule := x*nSolution;

        amountOfParticles := amountOfFreeMolecule/(1 - SelfClustering_K*x);

        specificAmountOfSubstance := amountOfParticles/mass;

        //specificAmountOfSubstance := 1/((SelfClustering_K + 1)*definition.MolarWeight);
      end if;

      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfParticles;

    function specificAmountOfFreeBaseMolecule
      "Amount of substance free base molecule per mass of the substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
        "Amount of substance free base molecule per substance mass";

    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K,amountOfBaseMolecules,x;
    algorithm
      if not selfClustering(definition) then
        specificAmountOfFreeBaseMolecule := 1/molarMassOfBaseMolecule(definition);
      else
        SelfClustering_dG :=selfClusteringBondEnthalpy(definition) - solution.T*
          selfClusteringBondEntropy(definition);

        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));

        amountOfBaseMolecules:=mass/molarMassOfBaseMolecule(definition);
        x:=((2*SelfClustering_K+nSolution/amountOfBaseMolecules) -
         sqrt((4*SelfClustering_K*nSolution/amountOfBaseMolecules)+
         (nSolution/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        specificAmountOfFreeBaseMolecule := (x*nSolution)/mass;

      end if;
      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfFreeBaseMolecule;

   function specificEnthalpy
     "Specific molar enthalpy of the substance with electric potential dependence"
     import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy
      "Specific enthalpy";

   protected
     Modelica.Units.SI.MolarEnergy SelfClustering_dG;
     Real SelfClustering_K;
   algorithm
     if selfClustering(definition) then
       SelfClustering_dG := selfClusteringBondEnthalpy(definition) - solution.T*
         selfClusteringBondEntropy(definition);
       SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));
     end if;

     specificEnthalpy := (molarEnthalpy(definition,solution) +
       (if selfClustering(definition) then
       selfClusteringBondEnthalpy(definition)*SelfClustering_K/(
       SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(definition);

   /*algorithm 

  specificEnthalpy := molarEnthalpy(
    definition,
    solution)/
    molarMassOfBaseMolecule(definition);*/
   end specificEnthalpy;

   function specificVolume "Specific volume of the substance"
     import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.SpecificVolume specificVolume "Specific volume";

   protected
     Modelica.Units.SI.MolarEnergy SelfClustering_dG;
     Real SelfClustering_K;
   algorithm
     if selfClustering(definition) then
       SelfClustering_dG := selfClusteringBondEnthalpy(definition) - solution.T*
         selfClusteringBondEntropy(definition);
       SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));
     end if;

     specificVolume := (molarVolume(definition,solution) +
     (if selfClustering(definition) then
       selfClusteringBondVolume(definition)*SelfClustering_K/(
       SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(definition);

   end specificVolume;

    function specificHeatCapacityCp
    "Specific heat capacity at constant pressure"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.SpecificHeatCapacity specificHeatCapacityCp
      "Specific heat capacity at constant pressure";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
    algorithm
      if selfClustering(definition) then
        SelfClustering_dG := selfClusteringBondEnthalpy(definition) - solution.T*
          selfClusteringBondEntropy(definition);
        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));
      end if;

      specificHeatCapacityCp := (molarHeatCapacityCp(
          definition,solution) + (if selfClustering(definition) then
        selfClusteringBondHeatCapacityCp(definition)*SelfClustering_K/(
        SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(definition);

        //TODO: + selfClusteringBondEnthalpy * der(K/(K + 1))/der(solution.T) .. if (selfClusteringBondHeatCapacityCp!=0)

    end specificHeatCapacityCp;

   function temperature
    "Temperature of the substance from its enthalpy"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
     input Modelica.Units.SI.MolarEnthalpy h "Molar enthalpy";

     input Modelica.Units.SI.Pressure p=100000 "Pressure";
     input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
     input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Temperature T "Temperature";

   protected
      function f_nonlinear "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
        extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
        input Interfaces.DataRecord data "Ideal gas data";
        input Modelica.Units.SI.MolarEnthalpy h "Molar enthalpy";
      algorithm
        y := H_T(data=data, T=u) - h;
      end f_nonlinear;
   algorithm
      T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function f_nonlinear(data=definition.data, h=h), 200, 6000);

   /*protected 
  Modelica.Units.SI.SpecificEnthalpy baseSpecificEnthalpy;
algorithm 

  baseSpecificEnthalpy := specificEnthalpy(
      definition,
      298.15,
      p,
      v,
      I);

  T := 298.15 + (h - baseSpecificEnthalpy)/specificHeatCapacityCp(
  substanceData); */
   end temperature;

   function solution_temperature
    "Temperature of the solution from specific enthalpy and mass fractions of substances"
     import Chemical;

       extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition[:] "Definition of substances";
    input Modelica.Units.SI.MolarEnthalpy h
      "Molar enthalpy of solution";
    input Modelica.Units.SI.MoleFraction x[:]
      "Mole fractions of substances";

    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Temperature T "Temperature";

   protected
     Chemical.Interfaces.Definition solutionDefinition=x*definition;

   algorithm
       T :=temperature(
           solutionDefinition,
           h,
           p,
           v,
           I);
   end solution_temperature;

   function density
        "Return density of the substance in the solution"
     import Chemical;
          extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.Density density "Density";

   algorithm
     density :=if definition.data.phase == Chemical.Interfaces.Phase.Gas then (definition.data.MM*solution.p)/(R*solution.T) else definition.data.MM/(definition.data.VmBase+definition.data.VmExcess);

   end density;

   function H_T "Molar enthalpy of the substance in electroneutral solution"
      extends Modelica.Icons.Function;
      import Modelica.Math;
      input Interfaces.DataRecord data "Data record of the substance";
      input Modelica.Units.SI.Temperature T "Temperature";

      output Modelica.Units.SI.MolarEnthalpy H "Molar enthalpy";
   algorithm
       H :=
       smooth(0,(if T < Tlimit then R*((-data.alow[1] + T*(data.blow[1] +
       data.alow[2]*Math.log(T) + T*(1.*data.alow[3] + T*(0.5*data.alow[4] +
       T*(1/3*data.alow[5] + T*(0.25*data.alow[6] + 0.2*data.alow[7]*T))))))
       /T) else R*((-data.ahigh[1] + T*(data.bhigh[1] + data.ahigh[2]*
       Math.log(T) + T*(1.*data.ahigh[3] + T*(0.5*data.ahigh[4] + T*(1/3*data.ahigh[5] +
       T*(0.25*data.ahigh[6] + 0.2*data.ahigh[7]*T))))))/T)) +
       data.H0);
       //(if   exclEnthForm then -data.Hf else 0.0) +
       //(if (refChoice == Choices.ReferenceEnthalpy.ZeroAt0K) then data.H0 else 0.0) +
       //(if refChoice == Choices.ReferenceEnthalpy.UserDefined then h_off else 0.0))

   end H_T;
    annotation (Documentation(revisions="<html>
<p><i>2015-2016</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end Properties;

  constant Modelica.Units.SI.Temperature Tlimit=1000 "Temperature limit between low and high data sets";

operator record DataRecord "Coefficient data record for chemical definitions based on NASA source"

  Modelica.Units.SI.MolarMass MM "Molar mass";
  Modelica.Units.SI.MolarEnthalpy Hf "Enthalpy of formation at 298.15K, 1bar";
  Modelica.Units.SI.MolarEnthalpy H0 "H0(298.15K, 1bar) - H0(0K, 1bar)";

  Real alow[7] "Low temperature coefficients a at 298.15K, 1bar";
  Real blow[3] "Low temperature constants b at 298.15K, 1bar";
  Real ahigh[7] "High temperature coefficients a at 298.15K, 1bar";
  Real bhigh[3] "High temperature constants b at 298.15K, 1bar";

  Modelica.Units.SI.ChargeNumberOfIon z
  "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

  Chemical.Interfaces.Phase phase "State of matter";

  // following is used only if phase is not Chemical.Interfaces.Phase.Gas:
  Modelica.Units.SI.MolarVolume VmBase "Base molar volume at 298.15K, 1bar (molar volume Vm = VmBase + VmExcess, density = MM/Vm)";
  Modelica.Units.SI.MolarVolume VmExcess "Excess of molar volume at 298.15K, 1bar (activity coeficient = exp(VmExcess/VmBase)";

   encapsulated operator 'constructor'
    import ModelicaDataRecord=Modelica.Media.IdealGases.Common.DataRecord;
    import DataRecord=Chemical.Interfaces.DataRecord;
    import PhaseType=Chemical.Interfaces.Phase;

      //    constant Real R=1.380649e-23*6.02214076e23;
      //    constant Real T0=298.15 "Base temperature";
      //    constant Real logT0 = 5.697596715569114904552663960891;
      //    constant Real Tlimit=1000;
      //    constant Real logTlimit=6.9077552789821370520539743640531;
      //    constant Real p0=100000 "Base pressure";

     function fromModelicaDataRecord
      input ModelicaDataRecord mdata "Mass based data record";
      input Real Gf=0 "Gibbs energy of formation at 298.15K, 1bar";
      input Real z=0 "Charge number";
      input PhaseType phase=PhaseType.Gas "State of matter";
      input Real Vm=if (phase == PhaseType.Gas) then (1.380649e-23*6.02214076e23)*298.15/100000 else 0.001*mdata.MM "Molar volume";
      input Real gamma=1 "Activity coefficient";

      output DataRecord result(
                  MM=mdata.MM,
                  Hf=mdata.Hf*mdata.MM,
                  H0=mdata.H0*mdata.MM,
                  alow=(mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23)) * mdata.alow,
                  blow={(mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23)) * mdata.blow[1],
                        (mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23)) * mdata.blow[2],
                         (((mdata.Hf*mdata.MM)-Gf)/298.15)/(1.380649e-23*6.02214076e23) - (mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23))*
                         (- 0.5*mdata.alow[1]/(298.15*298.15) - mdata.alow[2]/298.15 +
                         mdata.alow[3]*5.697596715569114904552663960891 + 298.15*(mdata.alow[4] +
                         298.15*(0.5*mdata.alow[5] + 298.15*(1/3*mdata.alow[6] +
                         0.25*mdata.alow[7]*298.15))))},
                        /*blow[3] = S(298.15)/(1.380649e-23*6.02214076e23) - rSlow(298.15), where S(T)=R* */
                  ahigh=(mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23)) * mdata.ahigh,

                  bhigh={(mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23)) * mdata.bhigh[1],
                         (mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23)) * mdata.bhigh[2],
                         ( ((((mdata.Hf*mdata.MM)-Gf)/298.15)/(1.380649e-23*6.02214076e23) - (mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23))*
                         (- 0.5*mdata.alow[1]/(298.15*298.15) - mdata.alow[2]/298.15 +
                         mdata.alow[3]*5.697596715569114904552663960891 + 298.15*(mdata.alow[4] +
                         298.15*(0.5*mdata.alow[5] + 298.15*(1/3*mdata.alow[6] +
                         0.25*mdata.alow[7]*298.15))))) +
                         (mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23))*(-0.5*mdata.alow[1]/(1000*1000)
                           - mdata.alow[2]/1000 + mdata.alow[3]*6.9077552789821370520539743640531
                           + 1000*(mdata.alow[4] + 1000*(0.5*mdata.alow[5]
                           + 1000*(1/3*mdata.alow[6] + 0.25*mdata.alow[7]*1000)))))
                         - (mdata.R_s*mdata.MM/(1.380649e-23*6.02214076e23))*
                         (- 0.5*mdata.ahigh[1]/(1000*1000) - mdata.ahigh[2]/1000 +
                         mdata.ahigh[3]*6.9077552789821370520539743640531 + 1000*(mdata.ahigh[4] +
                         1000*(0.5*mdata.ahigh[5] + 298.15*(1/3*mdata.ahigh[6] +
                         0.25*mdata.ahigh[7]*1000))))},
                         /*bhigh[3] = Slow(1000)/(1.380649e-23*6.02214076e23) - rShigh(1000) */


                  z=z,
                  phase=phase,
                  VmBase=Vm/(1+log(gamma)),
                  VmExcess=Vm*log(gamma)/(1+log(gamma)))
                           "Molar based data record";
     algorithm
       assert(mdata.Tlimit==1000, "Tlimit must be 1000K!");
      annotation (Inline=true);
     end fromModelicaDataRecord;


    function fromValues
      input Real MM;
      input Real Hf;
      input Real H0;
      input Real alow[7];
      input Real blow[3];
      input Real ahigh[7];
      input Real bhigh[3];
      input Real z=0 "Charge number";
      input PhaseType phase=PhaseType.Incompressible "State of matter";
      input Real VmBase=if (phase == PhaseType.Gas) then (1.380649e-23*6.02214076e23)*298.15/100000 else 0.001*MM "Base molar volume";
      input Real VmExcess=0 "Excess molar volume";

      output DataRecord result(
          MM = MM,
          Hf = Hf,
          H0 = H0,
          alow = alow,
          blow = blow,
          ahigh = ahigh,
          bhigh = bhigh,
          z = z,
          phase=phase,
          VmBase=VmBase,
          VmExcess=VmExcess);
    algorithm
      annotation (Inline=true);
    end fromValues;

   end 'constructor';
  annotation (Documentation(info="<html>
<p>
This data record contains the coefficients for the
ideal gas equations according to:
</p>
<blockquote>
  <p>McBride B.J., Zehe M.J., and Gordon S. (2002): <strong>NASA Glenn Coefficients
  for Calculating Thermodynamic Properties of Individual Species</strong>. NASA
  report TP-2002-211556</p>
</blockquote>
<p>
The equations have the following structure:
</p>
<div><img src=\"modelica://Modelica/Resources/Images/Media/IdealGases/Common/singleEquations.png\"></div>
<p>
The polynomials for h(T) and s0(T) are derived via integration from the one for cp(T)  and contain the integration constants b1, b2 that define the reference specific enthalpy and entropy. For entropy differences the reference pressure p0 is arbitrary, but not for absolute entropies. It is chosen as 1 standard atmosphere (101325 Pa).
</p>
<p>
For most gases, the region of validity is from 200 K to 6000 K.
The equations are split into two regions that are separated
by Tlimit (usually 1000 K). In both regions the gas is described
by the data above. The two branches are continuous and in most
gases also differentiable at Tlimit.
</p>
</html>"));
end DataRecord;



  connector DefinitionInput = input Chemical.Interfaces.Definition;
  connector DefinitionOutput = output Chemical.Interfaces.Definition;
  connector SubstanceStateInput = input SubstanceState;
  connector SubstanceStateOutput = output SubstanceState;
  connector SolutionStateInput = input SolutionState;
  connector SolutionStateOutput = output SolutionState;

  model Total "Summation of all extensible properties per substance"


    SolutionPort solution
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    parameter Boolean ElectricGround = true
    "Is the solution electric potential equal to zero during simulation (if not useElectricPort)?"
      annotation (HideResult=true);

    Modelica.Blocks.Interfaces.RealInput pressure "pressure"
      annotation (Placement(transformation(extent={{-120,58},{-80,100}})));
    Modelica.Blocks.Interfaces.RealInput temperature "temperature"
      annotation (Placement(transformation(extent={{-120,-22},{-80,20}})));
    Modelica.Blocks.Interfaces.RealOutput volume_der "derivation of volume"
      annotation (Placement(transformation(extent={{100,80},{120,100}})));
    Modelica.Blocks.Interfaces.RealOutput enthalpy_der "derivation of enthalpy"
      annotation (Placement(transformation(extent={{100,50},{120,70}})));
    Modelica.Blocks.Interfaces.RealOutput gibbsEnergy "Gibbs Energy of solution"
      annotation (Placement(transformation(extent={{100,20},{120,40}})));
    Modelica.Blocks.Interfaces.RealOutput charge
                                             "electric charge"
      annotation (Placement(transformation(extent={{100,-40},{120,-20}})));
    Modelica.Blocks.Interfaces.RealOutput volume
                                             "volume"
      annotation (Placement(transformation(extent={{100,-70},{120,-50}})));
    Modelica.Blocks.Interfaces.RealOutput mass
                                             "mass"
      annotation (Placement(transformation(extent={{100,-100},{120,-80}})));
  equation

    solution.p =pressure;
    solution.T =temperature;

    volume_der + solution.dV = 0;
    enthalpy_der + solution.dH = 0;

     //aliases
    solution.G = gibbsEnergy;
    solution.Q = charge;
    solution.V = volume;
    solution.m = mass;

    //electric current to solution must be represented with some substances (e.g. mass flow of electrones)
    if ElectricGround then
      //Solution connected to ground has zero voltage. However, electric current from the solution can varies.
        solution.v = 0;
    end if;
    if (not ElectricGround) then
      //Electrically isolated solution has not any electric current from/to the solution. However, electric potential can varies.
        solution.i = 0;
    end if;

    //Extensive properties of the solution:

    // The extensive quantities here have not the real physical flows.
    // They hack the Kirchhof's flow equation to be counted as the sum from all connected substances in the solution.

    //amount of substances
    solution.n + solution.nj = 0; //total amount of solution is the sum of amounts of each substance

    //mass of substances
    solution.m + solution.mj = 0; //total mass of solution is the sum masses of each substance

    //Gibs energy
    solution.G + solution.Gj = 0; //total free Gibbs energy of solution is the sum of free Gibbs energies of each substance

    //enthalpy
  //  solution.H + solution.Hj = 0;  //total free enthalpy of solution is the sum of enthalpies of each substance

    //ionic strength (mole fraction based)
    solution.I + solution.Ij = 0; //total ionic strength of solution is the ionic strengths of each substance

    //electric charge
    solution.Q + solution.Qj = 0; //total electric charge of solution is the sum of charges of each substance

    //volume
    solution.V + solution.Vj = 0; //total volume of solution is the sum of volumes of each substance

                                                                                                      annotation (
      Documentation(revisions="<html>
<p>2015-2018 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<h4>amountOfSubstances = &int; MolarFlows</h4>
<h4>mass = &int; massChanges</h4>
<h4>volume = &int; volumeChanges</h4>
<h4>freeEnthalpy = &int; EnthalpyChanges</h4>
<h4>freeEntropy = &int; EntropyChanges</h4>
<h4>freeGibbsEnergy = &int; GibbsEnergyChanges</h4>
<p>Integration of all substances together into one homogenous mixture - the solution.</p>
</html>"));
  end Total;

  partial model PartialSolution
    "Base chemical solution as homogenous mixture of the substances (only pressure and electric potential are not defined)"


    outer Modelica.Fluid.System system "System wide properties";

    parameter Boolean ElectricGround = true
    "Is electric potential equal to zero?"
      annotation (Evaluate=true, choices(checkBox=true), Dialog(group="Environment relationships"));

  Modelica.Units.SI.Temperature temperature "Temperature";

  Modelica.Units.SI.Pressure pressure "Pressure";

  Modelica.Units.SI.Volume volume "Current volume of the solution";

  Modelica.Units.SI.Mass mass(stateSelect=StateSelect.prefer)
    "Current mass of the solution";

    Total total(ElectricGround=ElectricGround)
      annotation (Placement(transformation(extent={{74,-96},{94,-76}})));

  protected
  Modelica.Units.SI.Energy gibbsEnergy
    "Gibbs energy of the solution relative to start of the simulation";

  Modelica.Units.SI.HeatFlowRate heatFromEnvironment
    "External heat flow rate";

  Modelica.Units.SI.ElectricCharge charge
    "Current electric charge of the solution";

  Modelica.Units.SI.HeatFlowRate enthalpy_der "derivative of enthalpy";
  Modelica.Units.SI.VolumeFlowRate volume_der "derivative of volume";

  equation

    heatFromEnvironment = enthalpy_der;

    //total inputs - thermodynamic state
    total.pressure = pressure;
    total.temperature = temperature;

    //total outputs = extensible properties
    enthalpy_der = total.enthalpy_der;
    volume_der = total.volume_der;
    gibbsEnergy = total.gibbsEnergy;
    charge = total.charge;
    volume = total.volume;
    mass = total.mass;

  end PartialSolution;

  partial model PartialSolutionWithHeatPort
    "Chemical solution as homogenous mixture of the substances"

    extends Interfaces.PartialSolution(temperature(start=temperature_start));

  parameter Modelica.Units.SI.Temperature temperature_start=system.T_ambient
    "Initial temperature of the solution"
    annotation (Dialog(group="Initialization"));

    parameter Boolean useThermalPort = false "Is thermal port pressent?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Boolean ConstantTemperature = true
    "Is temperature constant (if not useThermalPort)?"
       annotation (Evaluate=true, choices(checkBox=true), Dialog(enable=not useThermalPort, group="Environment relationships"));

    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort(T=temperature,
        Q_flow=heatFromEnvironment) if useThermalPort annotation (Placement(
          transformation(extent={{-70,-90},{-50,-70}}), iconTransformation(
            extent={{-62,-104},{-58,-100}})));

  initial equation
    temperature = temperature_start;
  equation

    //thermal
    if (not useThermalPort) and ConstantTemperature then
      //Ideal thermal exchange between environment and solution to reach constant temperature
      der(temperature) = 0;
    end if;
    if (not useThermalPort) and (not ConstantTemperature) then
      //Thermally isolated without any thermal exchange with environment
      heatFromEnvironment = 0;
    end if;

                                                                                                      annotation (
      Icon(coordinateSystem(
            preserveAspectRatio=false, initialScale=1, extent={{-100,-100},{
            100,100}})),
      Documentation(revisions="<html>
<p>2018 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<h4>amountOfSolution = &sum; amountOfSubstances</h4>
<h4>mass = &sum; massOfSubstances</h4>
<h4>volume = &sum; volumeOfSubstances</h4>
<h4>freeGibbsEnergy = &sum; freeGibbsEnergiesOfSubstances</h4>
<p>To calculate the sum of extensive substance's properties is misused the Modelica \"flow\" prefix even there are not real physical flows. </p>
</html>"));
  end PartialSolutionWithHeatPort;

  partial model ConditionalSubstanceFlow
  "Input of substance molar flow vs. parametric substance molar flow"

    parameter Boolean useSubstanceFlowInput = false
    "=true, if substance flow is provided via input"
    annotation(HideResult=true, choices(checkBox=true));

  parameter Modelica.Units.SI.MolarFlowRate SubstanceFlow=0
    "Volumetric flow of Substance if useSubstanceFlowInput=false"
    annotation (HideResult=true, Dialog(enable=not useSubstanceFlowInput));

    Modelica.Blocks.Interfaces.RealInput substanceFlow(start=SubstanceFlow, final unit="mol/s")=q
      if useSubstanceFlowInput
         annotation (HideResult=true,
         Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={40,40})));

  Modelica.Units.SI.MolarFlowRate q "Current Substance flow";
  equation
    if not useSubstanceFlowInput then
      q = SubstanceFlow;
    end if;

  end ConditionalSubstanceFlow;

  partial model PartialSolutionSensor

    import Chemical.Utilities.Types.SolutionChoice;

   // parameter SolutionChoice solutionFrom = Chemical.Utilities.Types.SolutionChoice.fromSubstrate "Chemical solution "
    parameter SolutionChoice solutionFrom = Chemical.Utilities.Types.SolutionChoice.FirstSubstrate "Chemical solution comes from?"
        annotation(HideResult=true, Dialog(group="Chemical solution (of products)"));

    parameter Chemical.Interfaces.SolutionState solutionParam = Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible) "Chemical solution state as Parameter"
      annotation (HideResult=true, Dialog(enable=(solutionFrom == SolutionChoice.Parameter), group="Chemical solution (of products)"));

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
          if (solutionFrom == SolutionChoice.SolutionPort) "To connect substance with solution, where is pressented"
      annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

      Chemical.Interfaces.SolutionState solutionState;

  protected

      Chemical.Interfaces.SolutionStateInput inputSubstrateSolution=solutionState if (solutionFrom == SolutionChoice.FirstSubstrate);

  equation

    if (solutionFrom == SolutionChoice.Parameter) then
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


     annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));

  end PartialSolutionSensor;

  partial model ConditionalKinetics
    "Input of kinetics coefficient vs. parametric kinetics coefficient"

    parameter Boolean useForwardRateInput=false
      "Forward rate coefficient from input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),
        Dialog(group="Conditional inputs", tab="Advanced"));

    parameter Real k_forward(unit="mol/s") = 1 "Forward rate coefficient (mole-fraction based)"
      annotation (HideResult=true, Dialog(enable=not useForwardRateInput));

    Modelica.Blocks.Interfaces.RealInput kfInput(start=k_forward, final unit="mol2.s-1.J-1")=
       kf if useForwardRateInput
       annotation ( HideResult=true, Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-60,40}),
                          iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-60,40})));

    Real kf(final unit="mol/s") "Current forward rate coefficient";

  equation
    if not useForwardRateInput then
      kf = k_forward;
    end if;

  end ConditionalKinetics;

  partial model SISO "Base Model with basic flow eqautions for SISO"

    import Chemical.Utilities.Types.InitializationMethods;

    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
      annotation(HideResult=true, Dialog(tab="Advanced"));

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the flow" annotation (HideResult=true, Dialog(tab="Advanced"));
    parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
      annotation(HideResult=true, Dialog(tab="Advanced"));
    parameter InitializationMethods initN_flow = Chemical.Utilities.Types.InitializationMethods.none "Initialization method for n_flow"
      annotation(HideResult=true, Dialog(tab= "Initialization"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_0=0 "Initial value for n_flow"
      annotation (HideResult=true, Dialog(tab="Initialization", enable=(initN_flow == InitializationMethods.state)));
    parameter Chemical.Utilities.Units.MolarFlowAcceleration n_acceleration_0=0 "Initial value for der(n_flow)"
      annotation (HideResult=true, Dialog(tab="Initialization", enable=(initN_flow == InitializationMethods.derivative)));


    Chemical.Interfaces.Rear rear( state_rearwards(u=u_rear_out, h=h_rear_out))
      annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

    Chemical.Interfaces.Fore fore( state_forwards(u=u_fore_out, h=h_fore_out))
      annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));


    Modelica.Units.SI.MolarFlowRate n_flow(stateSelect=n_flowStateSelect)=rear.n_flow;

  protected

    // delta potential computation
    Modelica.Units.SI.ChemicalPotential du_fore;
    Modelica.Units.SI.ChemicalPotential du_rear;

    Modelica.Units.SI.MoleFraction Px_rear, Sx_fore;
  /*  Modelica.Units.SI.Concentration c_in, c_out;
  Modelica.Units.SI.Molality b_in, b_out;
  Modelica.Units.SI.MassFraction X_in, X_out;
*/

    outer Chemical.DropOfCommons dropOfCommons;

    Chemical.Interfaces.SubstanceStateInput state_rear_in "Input substance state in forwards direction";
    Chemical.Interfaces.SubstanceStateInput state_fore_in "Input substance state in rearwards direction";

    //Chemical.Interfaces.SubstanceState state_out "Input substance state in rearwards direction";

    // input state quantities
    //Modelica.Units.SI.ChemicalPotential u_rear_in=rear.state_forwards.u "Chemical potential of substance entering";
    //Modelica.Units.SI.ChemicalPotential u_fore_in=fore.state_rearwards.u "Chemical potential of substance entering";
    //Modelica.Units.SI.MolarEnthalpy h_rear_in=rear.state_forwards.h "Enthalpy of substance enetering";
    //Modelica.Units.SI.MolarEnthalpy h_fore_in=fore.state_rearwards.h "Enthalpy of substance enetering";

    //outlet state quantities
    Modelica.Units.SI.ChemicalPotential u_rear_out "Chemical potential of substance exiting";
    Modelica.Units.SI.ChemicalPotential u_fore_out "Chemical potential of substance exiting";
    Modelica.Units.SI.MolarEnthalpy h_rear_out "Enthalpy of substance exiting";
    Modelica.Units.SI.MolarEnthalpy h_fore_out "Enthalpy of substance exiting";

    Modelica.Units.SI.ChemicalPotential uPure_substrate "Electro-chemical potential of pure substance entering";
    Modelica.Units.SI.ChemicalPotential uPure_product "Electro-chemical potential of pure substance exiting";

    //Chemical.Utilities.Units.URT duRT_fore, duRT_rear;

    Real Kx;

  initial equation
    if initN_flow == InitializationMethods.state then
      n_flow = n_flow_0;
    elseif initN_flow == InitializationMethods.derivative then
      der(n_flow) = n_acceleration_0;
    elseif initN_flow == InitializationMethods.steadyState then
      der(n_flow) = 0;
    end if;

  equation

    connect(rear.state_forwards, state_rear_in);
    connect(fore.state_rearwards, state_fore_in);

    du_fore = rear.state_forwards.u - fore.state_forwards.u;
    du_rear = fore.state_rearwards.u - rear.state_rearwards.u;

    /*
  du_fore = (s * substrates.state_forwards.u) - (p * products.state_forwards.u);
  du_rear = (p * products.state_rearwards.u) - (s * substrates.state_rearwards.u);
  */


    fore.state_forwards.h = rear.state_forwards.h;
    rear.state_rearwards.h = fore.state_rearwards.h;


  //  duRT_fore = (rear.state_forwards.u / (Modelica.Constants.R*rear.solution.T)) - (fore.state_forwards.u / (Modelica.Constants.R*fore.solution.T));
  //  duRT_rear = (rear.state_rearwards.u / (Modelica.Constants.R*rear.solution.T)) - (fore.state_rearwards.u / (Modelica.Constants.R*fore.solution.T));


    Sx_fore = exp(((rear.state_forwards.u - uPure_substrate)./(Modelica.Constants.R*rear.solution_forwards.T)));
    Px_rear = exp(((fore.state_rearwards.u - uPure_product)./(Modelica.Constants.R*fore.solution_rearwards.T)));
    Kx = exp(uPure_product/(Modelica.Constants.R*fore.solution_rearwards.T)-uPure_substrate/(Modelica.Constants.R*rear.solution_forwards.T));


    uPure_substrate = Chemical.Interfaces.Properties.electroChemicalPotentialPure(
      rear.definition,
      rear.solution_forwards);
    uPure_product = Chemical.Interfaces.Properties.electroChemicalPotentialPure(
      fore.definition,
      fore.solution_rearwards);



    fore.n_flow + rear.n_flow = 0;
    fore.r = rear.r - der(n_flow) * L;






    annotation (Documentation(info="<html>
<u>Interface class for all components with one fore and one rear port and a massflow without a mass storage between.</u>
<u>This class already implements the equations that are common for such components, namly the conservation of mass, the intertance equation, as well as the clipping of u_fore to u_min. </u>
<u>If u_fore should be lower the u_min, the remaining potential drop is added on the difference in inertial potential r, basically accelerating or decelerating the massflow. </u>
<u>The component offers different initialization methods for the massflow, as well as several parameters used in the equations above. </u>
<u>The clipping of the massflow can be turned off (this should be done by the modeler as a final modificator while extending to hide this option from the enduser).</u>
</html>"));
  end SISO;

end Interfaces;
