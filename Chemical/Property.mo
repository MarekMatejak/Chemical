within Chemical;
package Property "Calculation of chemical properties from definition and solution state"

  constant Real R=1.380649e-23*6.02214076e23;
  constant Real T0=298.15 "Base temperature";
  constant Real p0=100000 "Base pressure";

  model BaseProperties "Base properties of the substance"

    parameter Boolean FixedDefinition "definition==definitionParam";

    parameter Definition definitionParam "used only if FixedDefinition to help initialization";

    parameter Modelica.Units.SI.Mass m_start "Start value for mass of the substance";

    parameter Boolean SolutionObserverOnly = false "True if substance does not affect the solution";

    Interfaces.InputDefinition definition "Definition of the substance";

    Interfaces.InputSolutionState solutionState "State of the solution";

    InputAmountOfSubstance amountOfBaseMolecules
      "Amount of base molecules inside all clusters in compartment";

    InputMolarFlowRate n_flow "Molar change of the substance";

    InputHeatFlowRate h_flow "Substance enthaply change";

    Modelica.Units.SI.MoleFraction x "Mole fraction of the substance";

    Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
      "Molar concentration of particles";

    Modelica.Units.SI.MassConcentration M(displayUnit="mg/l")
          "Mass concentration";

    Modelica.Units.SI.Molality b(displayUnit="mmol/kg")
          "Molality";

    Modelica.Units.SI.MassFraction X "Mass fraction";

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

   // Local connector definition, used for equation balancing check

    connector InputMolarFlowRate = input Modelica.Units.SI.MolarFlowRate
      "Molar flow rate as input signal connector";
    connector InputHeatFlowRate = input Modelica.Units.SI.HeatFlowRate
      "Heat flow rate as input signal connector";

    connector InputSubstanceData = input Definition
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
                                     Chemical.Interfaces.SolutionStateParameters(
                                     T=system.T_ambient,
                                     p=system.p_ambient)))
      "Amount of free molecules not included inside any clusters in compartment";

    Modelica.Units.SI.AmountOfSubstance amountOfParticles(start=
         m_start*specificAmountOfParticles(
                                     definitionParam,
                                     Chemical.Interfaces.SolutionStateParameters(
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
   VmPure = molarVolumePure(definition,solutionState);
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

    x = amountOfFreeMolecule/solutionState.n "mole fraction [mol/mol]";

    c = amountOfParticles/solutionState.V "concentration [mol/m3]";

    M = amountOfParticles*molarMassOfBaseMolecule(definition)/solutionState.V "mass concentration [kg/m3]";

    b = amountOfParticles/solutionState.m "molality [mol/kg]";

    X  = amountOfParticles*molarMassOfBaseMolecule(definition)/solutionState.m "mass fraction [kg/kg]";

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
  end BaseProperties;

 function processData
   "Process changes of Gibbs energy, enthalpy, volume and heat capacity (products - reactants)"
      extends Modelica.Icons.Function;
  input Real K "Process dissociation constant (mole-fraction based) at 25°C,1bar";
  input Modelica.Units.SI.MolarEnergy dH=0 "Process molar enthalpy change at 25°C,1bar";
  input Modelica.Units.SI.MolarHeatCapacity dCp=0 "Process molar heat capacity change at 25°C,1bar";
  input Modelica.Units.SI.MolarVolume dVm=0 "Process molar volume change at 25°C,1bar";
  input Interfaces.PhaseType phase=Interfaces.PhaseType.Incompressible "State of matter";
  output Definition processData "Data record of process changes";
 algorithm
     processData :=  Definition(
       MM = 0,
       z = 0,
       DfG = -Modelica.Constants.R*T0*log(K),
       DfH = dH,
       Cp = dCp,
       Vm = dVm,
       phase = phase);
       //Name="processData",

 end processData;

 function firstProductDefinition
     "Return formation definition of the first product of chemical process"
      extends Modelica.Icons.Function;
  //input
  input Modelica.Units.SI.StoichiometricNumber s[:] "Stoichiometric reaction coefficient for substrates [nS]";
  input Modelica.Units.SI.StoichiometricNumber p[:] "Stoichiometric reaction coefficient for products [nP]";
  input Definition process "Data record of process changes";
  input Definition substrates[:] "Substrates definitions [nS]";
  input Definition products[:] "Other products definitions [nP-1]";
  output Definition firstProductDefinition "Definition of the first product in process";
 algorithm
      firstProductDefinition := (1/p[1]) * (p[2:end]*products - s*substrates - process);
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

 function activityCoefficient
  "Return activity coefficient of the substance in the solution"
    extends Modelica.Icons.Function;
    input Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
    output Real activityCoefficient "Activity Coefficient";
 algorithm
     activityCoefficient := 1;
 end activityCoefficient;

 function chargeNumberOfIon
  "Return charge number of the substance in the solution"
    extends Modelica.Icons.Function;
    input Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.ChargeNumberOfIon chargeNumberOfIon
    "Charge number of ion";
 algorithm
     chargeNumberOfIon := definition.data.z;
 end chargeNumberOfIon;

 function molarEnthalpyElectroneutral
  "Molar enthalpy of the substance in electroneutral solution"
    extends Modelica.Icons.Function;
    import Modelica.Math;
    input Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.MolarEnthalpy molarEnthalpyElectroneutral
    "Molar enthalpy";
 algorithm
     //Molar enthalpy:
     // - temperature and pressure shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - dp*Vm)
     //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - (dp/dT)*Vm. As a result dH = dT*Cv - dp*Vm for incompressible substances.

     molarEnthalpyElectroneutral :=
     smooth(0,(if solution.T < definition.data.Tlimit then R*((-definition.data.alow[1] + solution.T*(definition.data.blow[1] +
     definition.data.alow[2]*Math.log(solution.T) + solution.T*(1.*definition.data.alow[3] + solution.T*(0.5*definition.data.alow[4] +
     solution.T*(1/3*definition.data.alow[5] + solution.T*(0.25*definition.data.alow[6] + 0.2*definition.data.alow[7]*solution.T))))))
     /solution.T) else R*((-definition.data.ahigh[1] + solution.T*(definition.data.bhigh[1] + definition.data.ahigh[2]*
     Math.log(solution.T) + solution.T*(1.*definition.data.ahigh[3] + solution.T*(0.5*definition.data.ahigh[4] + solution.T*(1/3*definition.data.ahigh[5] +
     solution.T*(0.25*definition.data.ahigh[6] + 0.2*definition.data.ahigh[7]*solution.T))))))/solution.T)) +
     definition.data.H0);
     //(if   exclEnthForm then -definition.data.Hf else 0.0) +
     //(if (refChoice == Choices.ReferenceEnthalpy.ZeroAt0K) then definition.data.H0 else 0.0) +
     //(if refChoice == Choices.ReferenceEnthalpy.UserDefined then h_off else 0.0))

      /*substanceData.DfH + (solution.T - 298.15)*
      substanceData.Cp;*/
     //   - (p - 100000) * molarVolumePure(substanceData,solution.T,p,v,I);
 end molarEnthalpyElectroneutral;

 function molarEnthalpy
  "Molar enthalpy of the substance with electric potential dependence"
    extends Modelica.Icons.Function;
    input Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.MolarEnthalpy molarEnthalpy "Molar enthalpy";
 algorithm
    molarEnthalpy := molarEnthalpyElectroneutral(definition,solution) +
         Modelica.Constants.F*chargeNumberOfIon(definition,solution)*solution.v;
    annotation (Inline=true, smoothOrder=2);
 end molarEnthalpy;

 function molarEntropyPure
  "Molar entropy of the pure substance"
    extends Modelica.Icons.Function;
    import Modelica.Math;
    input Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";

  output Modelica.Units.SI.MolarEntropy molarEntropyPure
    "Molar entropy of the pure substance";

 algorithm
     //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
     //+ substanceData.Cv*log(solution.T/298.15);

     //Molar entropy shift:
     // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = solution.T*dS (small amount of added heat energy)
     // - pressure shift: with constant molar volume at constant temperature Vm*dP = -solution.T*dS (small amount of work)
     molarEntropyPure := if solution.T < definition.data.Tlimit then R*(definition.data.blow[2] - 0.5*definition.data.alow[
     1]/(solution.T*solution.T) - definition.data.alow[2]/solution.T + definition.data.alow[3]*Math.log(solution.T) + solution.T*(
     definition.data.alow[4] + solution.T*(0.5*definition.data.alow[5] + solution.T*(1/3*definition.data.alow[6] + 0.25*definition.data.alow[
     7]*solution.T)))) else R*(definition.data.bhigh[2] - 0.5*definition.data.ahigh[1]/(solution.T*solution.T) - definition.data.
     ahigh[2]/solution.T + definition.data.ahigh[3]*Math.log(solution.T) + solution.T*(definition.data.ahigh[4]
      + solution.T*(0.5*definition.data.ahigh[5] + solution.T*(1/3*definition.data.ahigh[6] + 0.25*definition.data.ahigh[7]*solution.T))))
     + (if definition.data.phase==Interfaces.PhaseType.Gas then R*log(solution.p/100000)
        else (molarVolumePure(definition,solution)/solution.T)*(solution.p - 100000));

     //For example at triple point of water should be solution.T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
     //As definition.data: http://www1.lsbu.ac.uk/water/water_phase_diagram.html
     //At solution.T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K

 end molarEntropyPure;

  function molarEntropy "Molar entropy of the substance in the solution"
        extends Modelica.Icons.Function;
    input Modelica.Units.SI.ChemicalPotential u
    "Electro-chemical potential of the substance";
    input Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.MolarEntropy molarEntropy "Molar entropy";
  algorithm
      molarEntropy :=  (u - molarEnthalpy(definition,solution))/solution.T;
  end molarEntropy;

 function chemicalPotentialPure "Chemical potential of the pure substance"
    extends Modelica.Icons.Function;
  input Definition definition "Definition of substance";
  input Interfaces.SolutionState solution "Chemical solution state";
  output Modelica.Units.SI.ChemicalPotential chemicalPotentialPure
    "Base chemical potential";
 algorithm
     chemicalPotentialPure :=  molarEnthalpyElectroneutral(definition,solution) - solution.T*molarEntropyPure(definition,solution);
 end chemicalPotentialPure;

 function electroChemicalPotentialPure
  "Electro-chemical potential of the pure substance"
    extends Modelica.Icons.Function;
  input Definition definition "Definition of substance";
  input Interfaces.SolutionState solution "Chemical solution state";
  output Modelica.Units.SI.ChemicalPotential
    electroChemicalPotentialPure "Base electro-chemical potential";
 algorithm
  electroChemicalPotentialPure := chemicalPotentialPure(
       definition,
       solution) + Modelica.Constants.F*chargeNumberOfIon(definition,solution)*solution.v;
 end electroChemicalPotentialPure;

 function molarVolumePure "Molar volume of the pure substance"
    extends Modelica.Icons.Function;
   input Definition definition "Definition of substance";
   input Interfaces.SolutionState solution "Chemical solution state";
  output Modelica.Units.SI.MolarVolume molarVolumePure "Molar volume";
 algorithm
   molarVolumePure := (if definition.data.phase==Interfaces.PhaseType.Gas then R*solution.T/solution.p else definition.data.Vm);
   //ideal gas
 end molarVolumePure;

 function molarVolumeExcess
  "Excess molar volume of the substance in the solution"
    extends Modelica.Icons.Function;
   input Definition definition "Definition of substance";
   input Interfaces.SolutionState solution "Chemical solution state";
  output Modelica.Units.SI.MolarVolume molarVolumeExcess
    "Excess molar volume of the substance in the solution";
 algorithm
    molarVolumeExcess := molarVolumePure(definition,solution)*
       log(activityCoefficient(definition,solution)); //zero if activityCoefficient==1
    annotation (Inline=true, smoothOrder=2);
 end molarVolumeExcess;

 function molarVolume "Molar volume of the substance"
    extends Modelica.Icons.Function;
   input Definition definition "Definition of substance";
   input Interfaces.SolutionState solution "Chemical solution state";

  output Modelica.Units.SI.MolarVolume molarVolume "Molar volume";
 algorithm
  molarVolume :=molarVolumePure(
       definition,
       solution) + molarVolumeExcess(
       definition,
       solution);
    annotation (Inline=true, smoothOrder=2);
 end molarVolume;

 function molarHeatCapacityCp
  "Molar heat capacity at constant pressure"
    extends Modelica.Icons.Function;
   input Definition definition "Definition of substance";
   input Interfaces.SolutionState solution "Chemical solution state";
  output Modelica.Units.SI.MolarHeatCapacity molarHeatCapacityCp
    "Molar heat capacity at constant pressure";
 algorithm
   molarHeatCapacityCp := smooth(0,if solution.T < definition.data.Tlimit then R*(1/(solution.T*solution.T)*(definition.data.alow[1] + solution.T*(
     definition.data.alow[2] + solution.T*(1.*definition.data.alow[3] + solution.T*(definition.data.alow[4] + solution.T*(definition.data.alow[5] + solution.T
     *(definition.data.alow[6] + definition.data.alow[7]*solution.T))))))) else R*(1/(solution.T*solution.T)*(definition.data.ahigh[1]
      + solution.T*(definition.data.ahigh[2] + solution.T*(1.*definition.data.ahigh[3] + solution.T*(definition.data.ahigh[4] + solution.T*(definition.data.
     ahigh[5] + solution.T*(definition.data.ahigh[6] + definition.data.ahigh[7]*solution.T))))))));
 end molarHeatCapacityCp;

 function molarMassOfBaseMolecule
    "Molar mass of base molecule of the substance"
    extends Modelica.Icons.Function;
    input Definition definition "Data record of substance";
  output Modelica.Units.SI.MolarMass molarMass "Molar mass";
 algorithm
   molarMass := definition.data.MM;
 end molarMassOfBaseMolecule;

 function selfClustering "returns true if substance molecules are joining together to clusters"
     extends Modelica.Icons.Function;
        input Definition definition "Data record of substance";
        output Boolean selfClustering;
 algorithm
   selfClustering:=false;
 end selfClustering;

 function selfClusteringBondEnthalpy
  "Enthalpy of joining two base molecules of the substance together to cluster"
     extends Modelica.Icons.Function;
        input Definition definition "Data record of substance";
  output Modelica.Units.SI.MolarEnthalpy selfClusteringEnthalpy;
 algorithm
   selfClusteringEnthalpy:=0;
 end selfClusteringBondEnthalpy;

 function selfClusteringBondEntropy
  "Entropy of joining two base molecules of the substance together to cluster"
     extends Modelica.Icons.Function;
        input Definition definition "Data record of substance";
  output Modelica.Units.SI.MolarEntropy selfClusteringEntropy;
 algorithm
   selfClusteringEntropy:=0;
 end selfClusteringBondEntropy;

 function selfClusteringBondVolume
     extends Modelica.Icons.Function;
        input Definition definition "Data record of substance";
  output Modelica.Units.SI.MolarVolume selfClusteringBondVolume;
 algorithm
   selfClusteringBondVolume:=0;
 end selfClusteringBondVolume;

 function selfClusteringBondHeatCapacityCp
    extends Modelica.Icons.Function;
        input Definition definition "Data record of substance";
  output Modelica.Units.SI.MolarHeatCapacity selfClusteringBondHeatCapacityCp;
 algorithm
   selfClusteringBondHeatCapacityCp:=0;
 end selfClusteringBondHeatCapacityCp;

  function specificAmountOfParticles
    "Amount of particles per mass of the substance"
    extends Modelica.Icons.Function;
    input Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
    input Modelica.Units.SI.Mass mass=1 "Mass of substance";
    input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
    output Real specificAmountOfSubstance(unit="mol/kg")
      "Amount of substance particles per its mass";
  algorithm
    specificAmountOfSubstance := 1/molarMassOfBaseMolecule(definition);
    annotation (Inline=true, smoothOrder=2);
  end specificAmountOfParticles;

  function specificAmountOfFreeBaseMolecule
    "Amount of substance free base molecule per mass of the substance"
    extends Modelica.Icons.Function;
    input Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
    input Modelica.Units.SI.Mass mass=1 "Mass of substance";
    input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
    output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
      "Amount of substance free base molecule per substance mass";
  algorithm
    specificAmountOfFreeBaseMolecule := 1/molarMassOfBaseMolecule(definition);
    annotation (Inline=true, smoothOrder=2);
  end specificAmountOfFreeBaseMolecule;

 function specificEnthalpy
   "Specific molar enthalpy of the substance with electric potential dependence"
    extends Modelica.Icons.Function;
  input Definition definition "Definition of substance";
  input Interfaces.SolutionState solution "Chemical solution state";

  output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy
    "Specific enthalpy";

 algorithm

   specificEnthalpy := molarEnthalpy(
     definition,
     solution)/
     molarMassOfBaseMolecule(definition);
 end specificEnthalpy;

 function specificVolume "Specific volume of the substance"
    extends Modelica.Icons.Function;
   input Definition definition "Definition of substance";
   input Interfaces.SolutionState solution "Chemical solution state";

  output Modelica.Units.SI.SpecificVolume specificVolume "Specific volume";

 algorithm

  specificVolume := molarVolume(
       definition,
       solution) /
     molarMassOfBaseMolecule(definition);
 end specificVolume;

  function specificHeatCapacityCp
  "Specific heat capacity at constant pressure"
    extends Modelica.Icons.Function;
    input Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
  output Modelica.Units.SI.SpecificHeatCapacity specificHeatCapacityCp
    "Specific heat capacity at constant pressure";

  algorithm

  specificHeatCapacityCp := molarHeatCapacityCp(
       definition,
       solution) /
     molarMassOfBaseMolecule(definition);
  end specificHeatCapacityCp;

 function temperature
  "Temperature of the substance from its enthalpy"
    extends Modelica.Icons.Function;
   input Definition definition "Definition of substance";
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
 end temperature;

 function solution_temperature
  "Temperature of the solution from specific enthalpy and mass fractions of substances"
     extends Modelica.Icons.Function;
    input Definition definition[:] "Definition of substances";
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
     Definition solutionDefinition=x*definition;

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
        extends Modelica.Icons.Function;
  input Definition definition "Definition of substance";
  input Interfaces.SolutionState solution "Chemical solution state";

  output Modelica.Units.SI.Density density "Density";

 algorithm
   density :=
     if definition.data.phase==Interfaces.PhaseType.Gas then
           (definition.data.MM*solution.p)/(R*solution.T)
     else definition.data.MM/definition.data.Vm;

 end density;

 function H_T "Molar enthalpy of the substance in electroneutral solution"
    extends Modelica.Icons.Function;
    import Modelica.Math;
    input Interfaces.DataRecord data "Data record of the substance";
    input Modelica.Units.SI.Temperature T "Temperature";

    output Modelica.Units.SI.MolarEnthalpy H "Molar enthalpy";
 algorithm
     H :=
     smooth(0,(if T < data.Tlimit then R*((-data.alow[1] + T*(data.blow[1] +
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
end Property;
