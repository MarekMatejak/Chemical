within Chemical;
package Interfaces "Chemical interfaces"
  extends Modelica.Icons.InterfacesPackage;

  partial package StateOfMatter "Abstract package for all state of matters"

   replaceable partial record SubstanceData
      "Definition data of the chemical substance"

   end SubstanceData;

   replaceable function activityCoefficient
    "Return activity coefficient of the substance in the solution"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

      output Real activityCoefficient "Activity Coefficient";
   end activityCoefficient;

   replaceable function chargeNumberOfIon
    "Return charge number of the substance in the solution"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.ChargeNumberOfIon chargeNumberOfIon
      "Charge number of ion";
   end chargeNumberOfIon;

   replaceable function molarEnthalpyElectroneutral
    "Molar enthalpy of the substance in electroneutral solution"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarEnthalpy molarEnthalpyElectroneutral
      "Molar enthalpy";
   end molarEnthalpyElectroneutral;

   function molarEnthalpy
    "Molar enthalpy of the substance with electric potential dependence"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarEnthalpy molarEnthalpy
      "Molar enthalpy";
   algorithm
      molarEnthalpy := molarEnthalpyElectroneutral(substanceData,T,p,v,I) +
           Modelica.Constants.F*chargeNumberOfIon(substanceData,T,p,v,I)*v;
      annotation (Inline=true, smoothOrder=2);
   end molarEnthalpy;

   replaceable function molarEntropyPure
    "Molar entropy of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarEntropy molarEntropyPure
      "Molar entropy of the pure substance";
   end molarEntropyPure;

    function molarEntropy "Molar entropy of the substance in the solution"
          extends Modelica.Icons.Function;
    input Modelica.Units.SI.ChemicalPotential u
      "Electro-chemical potential of the substance";
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarEntropy molarEntropy "Molar entropy";
    algorithm
        molarEntropy :=  (u - molarEnthalpy(substanceData,T,p,v,I))/T;
    end molarEntropy;

   function chemicalPotentialPure "Chemical potential of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.ChemicalPotential chemicalPotentialPure
      "Base chemical potential";
   algorithm
       chemicalPotentialPure :=  molarEnthalpyElectroneutral(substanceData,T,p,v,I) - T*molarEntropyPure(substanceData,T,p,v,I);
   end chemicalPotentialPure;

   function electroChemicalPotentialPure
    "Electro-chemical potential of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.ChemicalPotential
      electroChemicalPotentialPure "Base electro-chemical potential";
   algorithm
    electroChemicalPotentialPure := chemicalPotentialPure(
         substanceData,
         T,
         p,
         v,
         I) + Modelica.Constants.F*chargeNumberOfIon(substanceData,T,p,v,I)*v;
   end electroChemicalPotentialPure;

   replaceable function molarVolumePure "Molar volume of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.MolarVolume molarVolumePure "Molar volume";
   end molarVolumePure;

   function molarVolumeExcess
    "Excess molar volume of the substance in the solution"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.MolarVolume molarVolumeExcess
      "Excess molar volume of the substance in the solution";
   algorithm
      molarVolumeExcess := molarVolumePure(substanceData,T,p,v,I)*
         log(activityCoefficient(substanceData,T,p,v,I)); //zero if activityCoefficient==1
      annotation (Inline=true, smoothOrder=2);
   end molarVolumeExcess;

   replaceable function molarVolume "Molar volume of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarVolume molarVolume "Molar volume";
   algorithm
    molarVolume :=molarVolumePure(
         substanceData,
         T,
         p,
         v,
         I) + molarVolumeExcess(
         substanceData,
         T,
         p,
         v,
         I);
      annotation (Inline=true, smoothOrder=2);
   end molarVolume;

   replaceable function molarHeatCapacityCp
    "Molar heat capacity at constant pressure"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.MolarHeatCapacity molarHeatCapacityCp
      "Molar heat capacity at constant pressure";
   end molarHeatCapacityCp;

   replaceable function molarMassOfBaseMolecule
      "Molar mass of base molecule of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarMass molarMass "Molar mass";
   end molarMassOfBaseMolecule;

   replaceable function selfClustering "returns true if substance molecules are joining together to clusters"
       extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
          output Boolean selfClustering;
   algorithm
     selfClustering:=false;
   end selfClustering;

   replaceable function selfClusteringBondEnthalpy
    "Enthalpy of joining two base molecules of the substance together to cluster"
       extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarEnthalpy selfClusteringEnthalpy;
   algorithm
     selfClusteringEnthalpy:=0;
   end selfClusteringBondEnthalpy;

   replaceable function selfClusteringBondEntropy
    "Entropy of joining two base molecules of the substance together to cluster"
       extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarEntropy selfClusteringEntropy;
   algorithm
     selfClusteringEntropy:=0;
   end selfClusteringBondEntropy;

   replaceable function selfClusteringBondVolume
       extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarVolume selfClusteringBondVolume;
   algorithm
     selfClusteringBondVolume:=0;
   end selfClusteringBondVolume;

   replaceable function selfClusteringBondHeatCapacityCp
      extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarHeatCapacity selfClusteringBondHeatCapacityCp;
   algorithm
     selfClusteringBondHeatCapacityCp:=0;
   end selfClusteringBondHeatCapacityCp;

    replaceable function specificAmountOfParticles
      "Amount of particles per mass of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0 "Ionic strengh (mole fraction based)";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfSubstance(unit="mol/kg")
        "Amount of substance particles per its mass";
    algorithm
      specificAmountOfSubstance := 1/molarMassOfBaseMolecule(substanceData);
      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfParticles;

    replaceable function specificAmountOfFreeBaseMolecule
      "Amount of substance free base molecule per mass of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0 "Ionic strengh (mole fraction based)";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
        "Amount of substance free base molecule per substance mass";
    algorithm
      specificAmountOfFreeBaseMolecule := 1/molarMassOfBaseMolecule(substanceData);
      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfFreeBaseMolecule;

  /* replaceable function solution_temperature_
    "Temperature of the solution from specific enthalpy and mass fractions of substances"
     extends Modelica.Icons.Function;
    input SubstanceData substanceData[:] "Data record of substances";
  input Modelica.Units.SI.MolarEnthalpy h
    "Molar enthalpy of solution (x*substances_h)";
  input Modelica.Units.SI.MoleFraction x[:]
    "Mole fractions of substances";
  input Modelica.Units.SI.Pressure p=100000 "Pressure";
  input Modelica.Units.SI.ElectricPotential v=0
    "Electric potential of the substance";
  input Modelica.Units.SI.MoleFraction I=0
    "Ionic strengh (mole fraction based)";

  output Modelica.Units.SI.Temperature T "Temperature";
    annotation (__Dymola_DymolaStoredErrors(thetext="/*replaceable function solution_temperature_
  \"Temperature of the solution from specific enthalpy and mass fractions of substances\"
    extends Modelica.Icons.Function;
   input SubstanceData substanceData[:] \"Data record of substances\";
 input Modelica.Units.SI.MolarEnthalpy h
   \"Molar enthalpy of solution (x*substances_h)\";
 input Modelica.Units.SI.MoleFraction x[:]
   \"Mole fractions of substances\";
 input Modelica.Units.SI.Pressure p=100000 \"Pressure\";
 input Modelica.Units.SI.ElectricPotential v=0
   \"Electric potential of the substance\";
 input Modelica.Units.SI.MoleFraction I=0
   \"Ionic strengh (mole fraction based)\";

 output Modelica.Units.SI.Temperature T \"Temperature\";
"));
end solution_temperature_;
*/

   replaceable function specificEnthalpy
     "Specific molar enthalpy of the substance with electric potential dependence"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy
      "Specific enthalpy";

   algorithm

     specificEnthalpy := molarEnthalpy(
       substanceData,
       T,
       p,
       v,
       I)/
       molarMassOfBaseMolecule(substanceData);
   end specificEnthalpy;

   replaceable function specificVolume "Specific volume of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.SpecificVolume specificVolume "Specific volume";

   algorithm

    specificVolume := molarVolume(
         substanceData,
         T,
         p,
         v,
         I) /
       molarMassOfBaseMolecule(substanceData);
   end specificVolume;

    replaceable function specificHeatCapacityCp
    "Specific heat capacity at constant pressure"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.SpecificHeatCapacity specificHeatCapacityCp
      "Specific heat capacity at constant pressure";

    algorithm

    specificHeatCapacityCp := molarHeatCapacityCp(
         substanceData,
         T,
         p,
         v,
         I) /
       molarMassOfBaseMolecule(substanceData);
    end specificHeatCapacityCp;

   replaceable function temperature
    "Temperature of the substance from its enthalpy"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.SpecificEnthalpy h "Specific enthalpy";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Temperature T "Temperature";
   end temperature;

   replaceable function solution_temperature
    "Temperature of the solution from specific enthalpy and mass fractions of substances"
       extends Modelica.Icons.Function;
      input SubstanceData substanceData[:] "Data record of substances";
    input Modelica.Units.SI.SpecificEnthalpy h
      "Specific enthalpy of solution";
    input Modelica.Units.SI.MassFraction X[:]
      "Mass fractions of substances";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Temperature T "Temperature";
   end solution_temperature;

   replaceable function density
        "Return density of the substance in the solution"
          extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Density density "Density";
   end density;
    annotation (Documentation(revisions="<html>
<p><i>2015-2016</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end StateOfMatter;

  package Incompressible "Incompressible as basic state of matter"
    extends StateOfMatter;

    redeclare record extends SubstanceData "Base substance data"

      parameter Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa")=
           0.01801528 "Molar weight of the substance";

      parameter Modelica.Units.SI.ChargeNumberOfIon z=0
        "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

      parameter Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol")=
        DfG_25degC_1bar
        "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      parameter Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol")=
        DfH_25degC
        "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      parameter Modelica.Units.SI.ActivityCoefficient gamma=1
        "Activity coefficient of the substance";

      parameter Modelica.Units.SI.MolarHeatCapacity Cp=0
        "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";
      parameter String References[1]={""}
        "References of these thermodynamical values";

      parameter Modelica.Units.SI.MolarEnergy DfG_25degC_1bar(displayUnit="kJ/mol")=
           0 "Obsolete parameter use DfH instead"
        annotation (Dialog(tab="Obsolete"));

      parameter Modelica.Units.SI.MolarEnergy DfH_25degC(displayUnit="kJ/mol")=
           0 "Obsolete parameter use DfG instead"
        annotation (Dialog(tab="Obsolete"));

      parameter Boolean SelfClustering=false
        "Pure substance is making clusters (weak bonds between molecules)";

      parameter Modelica.Units.SI.ChemicalPotential SelfClustering_dH=0
        "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";
      //-20000
      parameter Modelica.Units.SI.MolarEntropy SelfClustering_dS=0
        "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";

      parameter Modelica.Units.SI.Density density(displayUnit="kg/dm3") = 1000
        "Density of the pure substance (default density of water at 25degC)";

      //      parameter Modelica.SIunits.MolarHeatCapacity Cv = Cp
      //      "Molar heat capacity of the substance at constant volume";

      annotation (preferredView="info", Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceData;

    redeclare function extends activityCoefficient
      "Return activity coefficient of the substance in the solution"
    algorithm
      activityCoefficient := substanceData.gamma;
    end activityCoefficient;

    redeclare function extends chargeNumberOfIon
      "Return charge number of the substance in the solution"
    algorithm
      chargeNumberOfIon := substanceData.z;
    end chargeNumberOfIon;

    redeclare function extends molarEnthalpyElectroneutral
      "Molar enthalpy of the pure electroneutral substance"
    algorithm
      //Molar enthalpy:
      // - temperature and pressure shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - dp*Vm)
      //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - (dp/dT)*Vm. As a result dH = dT*Cv - dp*Vm for incompressible substances.

      molarEnthalpyElectroneutral := substanceData.DfH + (T - 298.15)*
        substanceData.Cp;
      //   - (p - 100000) * molarVolumePure(substanceData,T,p,v,I);
    end molarEnthalpyElectroneutral;

    redeclare function extends molarEntropyPure
      "Molar entropy of the pure substance"
    algorithm
      //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
      //+ substanceData.Cv*log(T/298.15);

      //Molar entropy shift:
      // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
      // - pressure shift: with constant molar volume at constant temperature Vm*dP = -T*dS (small amount of work)
      molarEntropyPure := substanceData.Cp*log(T/298.15) - (molarVolumePure(
          substanceData,
          T,
          p,
          v,
          I)/T)*(p - 100000) + ((substanceData.DfH - substanceData.DfG)/298.15);

      //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
      //As data: http://www1.lsbu.ac.uk/water/water_phase_diagram.html
      //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
    end molarEntropyPure;

    redeclare function molarVolumePure
      "Molar volume of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.MolarVolume molarVolumePure "Molar volume";
    algorithm
      molarVolumePure := substanceData.MolarWeight/substanceData.density;
      //incompressible
    end molarVolumePure;

    redeclare function extends molarHeatCapacityCp
      "Molar heat capacity of the substance at constant pressure"
    algorithm
      molarHeatCapacityCp := substanceData.Cp;
    end molarHeatCapacityCp;

    redeclare function extends molarMassOfBaseMolecule
      "Molar mass of the substance"
    algorithm
      molarMass := substanceData.MolarWeight;
    end molarMassOfBaseMolecule;

    redeclare function selfClustering
      "returns true if substance molecules are joining together to clusters"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      output Boolean selfClustering;
    algorithm
      selfClustering := substanceData.SelfClustering;
    end selfClustering;

    redeclare function selfClusteringBondEnthalpy
      "Enthalpy of joining two base molecules of the substance together to cluster"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      output Modelica.Units.SI.MolarEnthalpy selfClusteringEnthalpy;
    algorithm
      selfClusteringEnthalpy := substanceData.SelfClustering_dH;
    end selfClusteringBondEnthalpy;

    redeclare function selfClusteringBondEntropy
      "Entropy of joining two base molecules of the substance together to cluster"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      output Modelica.Units.SI.MolarEntropy selfClusteringEntropy;
    algorithm
      selfClusteringEntropy := substanceData.SelfClustering_dS;
    end selfClusteringBondEntropy;

    redeclare replaceable function specificAmountOfParticles
    "Amount of substance particles per its mass"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfSubstance(unit="mol/kg") "Amount of substance particles per its mass";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
      Real amountOfBaseMolecules,amountOfFreeMolecule,amountOfParticles;
      Real x;
    algorithm
      if not selfClustering(substanceData) then
        specificAmountOfSubstance := 1/substanceData.MolarWeight;
      else
        SelfClustering_dG :=selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);

        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));

        amountOfBaseMolecules:=mass/substanceData.MolarWeight;
        x:=((2*SelfClustering_K+nSolution/amountOfBaseMolecules) -
         sqrt((4*SelfClustering_K*nSolution/amountOfBaseMolecules)+
         (nSolution/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        amountOfFreeMolecule := x*nSolution;

        amountOfParticles := amountOfFreeMolecule/(1 - SelfClustering_K*x);

        specificAmountOfSubstance := amountOfParticles/mass;

        //specificAmountOfSubstance := 1/((SelfClustering_K + 1)*substanceData.MolarWeight);
      end if;
    end specificAmountOfParticles;

    redeclare function specificAmountOfFreeBaseMolecule
      "Amount of substance free base molecule per mass of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0 "Ionic strengh (mole fraction based)";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance in solution";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
        "Amount of substance free base molecule per substance mass";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K,amountOfBaseMolecules,x;
    algorithm
      if not selfClustering(substanceData) then
        specificAmountOfFreeBaseMolecule := 1/substanceData.MolarWeight;
      else
        SelfClustering_dG :=selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);

        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));

        amountOfBaseMolecules:=mass/substanceData.MolarWeight;
        x:=((2*SelfClustering_K+nSolution/amountOfBaseMolecules) -
         sqrt((4*SelfClustering_K*nSolution/amountOfBaseMolecules)+
         (nSolution/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        specificAmountOfFreeBaseMolecule := (x*nSolution)/mass;

      end if;
      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfFreeBaseMolecule;

    redeclare replaceable function specificEnthalpy
      "Specific molar enthalpy of the substance with electric potential dependence"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy
        "Specific enthalpy";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
    algorithm
      if selfClustering(substanceData) then
        SelfClustering_dG := selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);
        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
      end if;

      specificEnthalpy := (molarEnthalpy(
          substanceData,
          T,
          p,
          v,
          I) + (if selfClustering(substanceData) then
        selfClusteringBondEnthalpy(substanceData)*SelfClustering_K/(
        SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(substanceData);

      annotation (Inline=true, smoothOrder=2);
    end specificEnthalpy;

    redeclare replaceable function specificVolume
      "Specific volume of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.SpecificVolume specificVolume
        "Specific volume";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
    algorithm
      if selfClustering(substanceData) then
        SelfClustering_dG := selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);
        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
      end if;

      specificVolume := (molarVolume(
          substanceData,
          T,
          p,
          v,
          I) + (if selfClustering(substanceData) then
        selfClusteringBondVolume(substanceData)*SelfClustering_K/(
        SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(substanceData);
    end specificVolume;

    redeclare replaceable function specificHeatCapacityCp
      "Specific heat capacity at constant pressure"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.SpecificHeatCapacity specificHeatCapacityCp
        "Specific heat capacity at constant pressure";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
    algorithm
      if selfClustering(substanceData) then
        SelfClustering_dG := selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);
        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
      end if;

      specificHeatCapacityCp := (molarHeatCapacityCp(
          substanceData,
          T,
          p,
          v,
          I) + (if selfClustering(substanceData) then
        selfClusteringBondHeatCapacityCp(substanceData)*SelfClustering_K/(
        SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(substanceData);

        //TODO: + selfClusteringBondEnthalpy * der(K/(K + 1))/der(T) .. if (selfClusteringBondHeatCapacityCp!=0)
    end specificHeatCapacityCp;

    redeclare function extends temperature
      "Temperature of substance from its enthalpy"
    protected
      Modelica.Units.SI.SpecificEnthalpy baseSpecificEnthalpy;
    algorithm

      baseSpecificEnthalpy := specificEnthalpy(
          substanceData,
          298.15,
          p,
          v,
          I);

      T := 298.15 + (h - baseSpecificEnthalpy)/specificHeatCapacityCp(
        substanceData);
    end temperature;

    redeclare function extends solution_temperature
      "Temperature of the solution from enthalpies os substances"
      // Modelica.Units.SI.MoleFraction x[size(X, 1)];
    protected
      Modelica.Units.SI.SpecificEnthalpy solution_h_base;
    /*  Modelica.Units.SI.SpecificHeatCapacity solution_Cp=sum(X[i]*
      substanceData[i].Cp/molarMassOfBaseMolecule(substanceData[i]) for
      i in 1:size(X, 1));*/
    algorithm
      solution_h_base := X*specificEnthalpy(
          substanceData,
          298.15,
          p,
          v,
          I);
      T := 298.15 + (h - solution_h_base)/(X*specificHeatCapacityCp(substanceData));
    end solution_temperature;

     redeclare function extends density
      "Return density of the substance in the solution"
     algorithm
      density := substanceData.density;
     end density;

    annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end Incompressible;

  package IdealGas "Ideal gas with constant heat capacity"
     extends StateOfMatter;

     redeclare record extends SubstanceData "Base substance data"

      parameter Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa")=
         0.01801528 "Molar weight of the substance";

      parameter Modelica.Units.SI.ChargeNumberOfIon z=0
      "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

      parameter Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol")=
         DfG_25degC_1bar
      "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      parameter Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol")=
         DfH_25degC
      "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      parameter Modelica.Units.SI.ActivityCoefficient gamma=1
      "Activity coefficient of the substance";

      parameter Modelica.Units.SI.MolarHeatCapacity Cp=0
      "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";
     parameter String References[1]={""}
       "References of these thermodynamical values";

      parameter Modelica.Units.SI.MolarEnergy DfG_25degC_1bar(displayUnit=
         "kJ/mol") = 0 "Obsolete parameter use DfH instead"
      annotation (Dialog(tab="Obsolete"));

      parameter Modelica.Units.SI.MolarEnergy DfH_25degC(displayUnit=
          "kJ/mol") = 0 "Obsolete parameter use DfG instead"
      annotation (Dialog(tab="Obsolete"));

     parameter Boolean SelfClustering = false "Pure substance is making clusters (weak bonds between molecules)";

      parameter Modelica.Units.SI.ChemicalPotential SelfClustering_dH=0
      "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";                                                                    //-20000
      parameter Modelica.Units.SI.MolarEntropy SelfClustering_dS=0
      "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";

      annotation ( preferredView = "info", Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
     end SubstanceData;

   redeclare function extends activityCoefficient
    "Return activity coefficient of the substance in the solution"
   algorithm
       activityCoefficient := substanceData.gamma;
   end activityCoefficient;

   redeclare function extends chargeNumberOfIon
    "Return charge number of the substance in the solution"
   algorithm
      chargeNumberOfIon := substanceData.z;
   end chargeNumberOfIon;

   redeclare function extends molarEnthalpyElectroneutral
    "Molar enthalpy of the pure substance in electroneutral solution"
   algorithm
       //Molar enthalpy:
       // - temperature shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - R*dT)
       //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - R. As a result dH = dT*(Cv+R) for ideal gas.
       //   And the relation with molar heat capacity at constant pressure as Cp=Cv+R makes dH = dT*Cp.
       molarEnthalpyElectroneutral := substanceData.DfH
         +(T-298.15)*(substanceData.Cp);
   end molarEnthalpyElectroneutral;

   redeclare function extends molarEntropyPure
    "Molar entropy of the pure substance"
   algorithm
     //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
     //+ (substanceData.Cp+Modelica.Constants.R)*log(T/298.15);

       //Molar entropy:
       // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
       // - pressure shift: to reach the ideal gas equation at constant temperature Vm*dP = -T*dS (small amount of work)
       molarEntropyPure := (substanceData.Cp)*log(T/298.15) - Modelica.Constants.R*log(p/100000) + ((substanceData.DfH
        - substanceData.DfG)/298.15);

       //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
       //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
   end molarEntropyPure;

   redeclare function extends molarVolumePure
    "Molar volume of the pure substance"
   algorithm
       molarVolumePure := Modelica.Constants.R*T/p; //ideal gas
   end molarVolumePure;

   redeclare function extends molarHeatCapacityCp
    "Molar heat capacity of the substance at constant pressure"
   algorithm
       molarHeatCapacityCp := substanceData.Cp;
   end molarHeatCapacityCp;

   redeclare function extends molarMassOfBaseMolecule "Molar mass of the substance"
   algorithm
       molarMass := substanceData.MolarWeight;
   end molarMassOfBaseMolecule;

   redeclare function extends temperature "Temperature of substance from its enthalpy"
   algorithm
        T := 298.15 + (h-specificEnthalpy(substanceData,298.15,p,v,I))/specificHeatCapacityCp(substanceData);
   end temperature;

   redeclare function extends solution_temperature
    "Temperature of the solution from enthalpies os substances"
   algorithm
        T := 298.15 + (h-X*specificEnthalpy(
             substanceData,
             298.15,
             p,
             v,
             I))/(X*specificHeatCapacityCp(substanceData));
   end solution_temperature;

   redeclare function extends density
        "Return density of the substance in the solution"
   algorithm
           density := substanceData.MolarWeight * p / (Modelica.Constants.R * T);
   end density;

    annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end IdealGas;

  package IdealGasMSL "Ideal gas from Modelica Standard Library 3.2"
    extends StateOfMatter;

    redeclare record SubstanceData

      parameter Modelica.Media.IdealGases.Common.DataRecord data=Modelica.Media.IdealGases.Common.SingleGasesData.N2 "Definition of the substance";

    parameter Modelica.Units.SI.ChargeNumberOfIon z=0
      "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

    end SubstanceData;

    redeclare function extends activityCoefficient
      "Return activity coefficient of the substance in the solution"
    algorithm
      activityCoefficient := 1;
      annotation (Inline=true, smoothOrder=2);
    end activityCoefficient;

    redeclare function extends chargeNumberOfIon
      "Return charge number of the substance in the solution"
    algorithm
      chargeNumberOfIon := substanceData.z;
      annotation (Inline=true, smoothOrder=2);
    end chargeNumberOfIon;

    redeclare function extends molarEnthalpyElectroneutral
      "Molar enthalpy of the pure substance in electroneutral solution"
    algorithm
      //Molar enthalpy:
      // - temperature shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - R*dT)
      //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - R. As a result dH = dT*(Cv+R) for ideal gas.
      //   And the relation with molar heat capacity at constant pressure as Cp=Cv+R makes dH = dT*Cp.
      molarEnthalpyElectroneutral := substanceData.data.MM*
        Modelica.Media.IdealGases.Common.Functions.h_T(
          substanceData.data,
          T,
          false,
          Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C);
      annotation (Inline=true, smoothOrder=2);
    end molarEnthalpyElectroneutral;

    redeclare function extends molarEntropyPure
      "Molar entropy of the pure substance"
    algorithm
      //molarEntropyPure := ((substanceData.data.DfH - substanceData.data.DfG)/298.15)
      //+ (substanceData.data.Cp+Modelica.Constants.R)*log(T/298.15);

      //Molar entropy:
      // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
      // - pressure shift: to reach the ideal gas equation at constant temperature Vm*dP = -T*dS (small amount of work)

      molarEntropyPure := substanceData.data.MM*(
        Modelica.Media.IdealGases.Common.Functions.s0_T(substanceData.data, T) -
        substanceData.data.R_s*log(p/100000));

      //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
      //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
      annotation (Inline=true, smoothOrder=2);
    end molarEntropyPure;

    redeclare function extends molarVolumePure
      "Molar volume of the pure substance"
    algorithm
      molarVolumePure := substanceData.data.MM*substanceData.data.R_s*T/p;
      //ideal gas
      annotation (Inline=true, smoothOrder=2);
    end molarVolumePure;

    redeclare function extends molarHeatCapacityCp
      "Molar heat capacity of the substance at constant pressure"
    algorithm
      molarHeatCapacityCp := substanceData.data.MM*
        Modelica.Media.IdealGases.Common.Functions.cp_T(substanceData.data, T);
      annotation (Inline=true, smoothOrder=2);
    end molarHeatCapacityCp;

    redeclare function extends molarMassOfBaseMolecule "Molar mass of the substance"
    algorithm
      molarMass := substanceData.data.MM;
      annotation (Inline=true, smoothOrder=2);
    end molarMassOfBaseMolecule;

    redeclare function extends temperature "Temperature of substance from its enthalpy"
    protected
       function f_nonlinear "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
         extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
         input Modelica.Media.IdealGases.Common.DataRecord data "Ideal gas data";
         input Modelica.Units.SI.SpecificEnthalpy h "Specific enthalpy";
       algorithm
         y := Modelica.Media.IdealGases.Common.Functions.h_T(data=data, T=u,
         exclEnthForm=false,refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C)
                - h;
       end f_nonlinear;

    algorithm
       T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
         function f_nonlinear(data=substanceData.data, h=h), 200, 6000);
    end temperature;

    redeclare function extends solution_temperature
    "Temperature of the solution from enthalpies os substances"
     // Modelica.Units.SI.MolarMass MM=x*substanceData.data.MM
     //   "molar mass of solution";
     // Modelica.Units.SI.MassFraction x_mass[:]=(x .* substanceData.data.MM) ./
     //     MM "mass fractions";
    protected
        Modelica.Media.IdealGases.Common.DataRecord solutionData=
           Modelica.Media.IdealGases.Common.DataRecord(
               name="solution_temperature",
               MM= 1/sum(X./substanceData.data.MM),
               Hf= X*substanceData.data.Hf,
               H0= X*substanceData.data.H0,
               Tlimit = X*substanceData.data.Tlimit,
               alow = X*substanceData.data.alow,
               blow = X*substanceData.data.blow,
               ahigh = X*substanceData.data.ahigh,
               bhigh = X*substanceData.data.bhigh,
               R_s = X*substanceData.data.R_s);
                     //),
            //sum through moles, not masses

    algorithm
        T := temperature(SubstanceData(data=solutionData,z=X*(substanceData.z./molarMassOfBaseMolecule(substanceData))),h,p,v,I);
    end solution_temperature;

    redeclare function extends density
      "Return density of the substance in the solution"
    algorithm
      density := p/(substanceData.data.R_s*T);
      annotation (Inline=true, smoothOrder=2);
    end density;
    annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end IdealGasMSL;

  package IdealGasShomate "Ideal gas based on Shomate equations"
     extends StateOfMatter;

   redeclare record extends SubstanceData
    "Base substance data based on Shomate equations http://old.vscht.cz/fch/cz/pomucky/fchab/Shomate.html"

    parameter Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa")=
         0.01801528 "Molar weight of the substance";

    parameter Modelica.Units.SI.ChargeNumberOfIon z=0
      "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

    parameter Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol")=
         DfG_25degC_1bar
      "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

    parameter Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol")=
         DfH_25degC
      "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

    parameter Modelica.Units.SI.ActivityCoefficient gamma=1
      "Activity coefficient of the substance";

    parameter Modelica.Units.SI.MolarHeatCapacity Cp=cp_25degC
      "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";
     parameter String References[1]={""}
       "References of these thermodynamical values";

    parameter Modelica.Units.SI.MolarEnergy DfG_25degC_1bar(displayUnit=
         "kJ/mol") = 0 "Obsolete parameter use DfH instead"
      annotation (Dialog(tab="Obsolete"));

    parameter Modelica.Units.SI.MolarEnergy DfH_25degC(displayUnit=
          "kJ/mol") = 0 "Obsolete parameter use DfG instead"
      annotation (Dialog(tab="Obsolete"));

     parameter Boolean SelfClustering = false "Pure substance is making clusters (weak bonds between molecules)";

    parameter Modelica.Units.SI.ChemicalPotential SelfClustering_dH=0
      "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";                                                                    //-20000
    parameter Modelica.Units.SI.MolarEntropy SelfClustering_dS=0
      "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";

        parameter Real B(unit="J.mol-1")=0 "Shomate parameter B";
        parameter Real C(unit="J.mol-1")=0 "Shomate parameter C";
        parameter Real D(unit="J.K.mol-1")=0 "Shomate parameter D";
        parameter Real E(unit="J.K2.mol-1")=0 "Shomate parameter E";
        parameter Real X=0 "Shomate parameter X";
        parameter Real A_(unit="J.K-1.mol-1")=0 "Shomate parameter A'";
        parameter Real E_(unit="K")=1e-8 "Shomate parameter E'";

        parameter Real cp_25degC(unit="J.K-1.mol-1") = 33.6
         "Obsolete parameter use Cp instead"
         annotation (Dialog(tab="Obsolete"));

      annotation (preferredView = "info", Documentation(revisions="<html>
<p><i>2016-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
   end SubstanceData;

   redeclare function extends activityCoefficient
    "Return activity coefficient of the substance in the solution"
   algorithm
       activityCoefficient := substanceData.gamma;
   end activityCoefficient;

   redeclare function extends chargeNumberOfIon
    "Return charge number of the substance in the solution"
   algorithm
      chargeNumberOfIon := substanceData.z;
   end chargeNumberOfIon;

   redeclare function extends molarEnthalpyElectroneutral
    "Molar enthalpy of the pure substance in electroneutral solution, where der(Hm)=cp*der(T)"
    protected
     parameter Real T0=298.15;
     Real t=T/1000;
     parameter Real A=substanceData.Cp
       - ((10^6 * substanceData.A_* exp(1000*substanceData.E_)/T0)) / ((-1 + exp((1000*substanceData.E_)/T0))^2 * T0^2)
       - (10^6 * substanceData.E)/T0^2 - 0.001*substanceData.B*T0 - 10^(-6) * substanceData.C * T0^2
       - 10^(-9) * substanceData.D * T0^3 - sqrt(1/1000)* T0^0.5 * substanceData.X;

     parameter Real H=substanceData.DfH
       - 1000*(substanceData.A_/((-1 + exp((1000*substanceData.E_)/T0))*substanceData.E_)
       - (1000*substanceData.E)/T0 + 0.001*A*T0
       + 5.*10^(-7)*substanceData.B*T0^2 + (1/3)*10^(-9)*substanceData.C*T0^3
       + 2.5*10^(-13)*substanceData.D*T0^4 + (1/1000)^(1.5)/1.5 * T0^1.5 * substanceData.X);

   algorithm
       //Molar enthalpy:
       // - temperature shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - R*dT)
       //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - R. As a result dH = dT*(Cv+R) for ideal gas.
       //   And the relation with molar heat capacity at constant pressure as Cp=Cv+R makes dH = dT*Cp.
       molarEnthalpyElectroneutral :=
       H + 1000*(A*t + substanceData.B*t^2/2 + substanceData.C*t^3/3
       + substanceData.D*t^4/4 - substanceData.E/t + substanceData.X*t^1.5/1.5
       + substanceData.A_/substanceData.E_/(exp(substanceData.E_/t) - 1));

   end molarEnthalpyElectroneutral;

   redeclare function extends molarEntropyPure
    "Molar entropy of the pure substance, where der(Sm) = cp*der(T)/T"
    protected
     parameter Real T0=298.15;
     Real t=T/1000;
     parameter Real A= substanceData.Cp
       - ((10^6 * substanceData.A_* exp(1000*substanceData.E_)/T0)) / ((-1 + exp((1000*substanceData.E_)/T0))^2 * T0^2)
       - (10^6 * substanceData.E)/T0^2 - 0.001*substanceData.B*T0 - 10^(-6) * substanceData.C * T0^2
       - 10^(-9) * substanceData.D * T0^3 - sqrt(1/1000)* T0^0.5 * substanceData.X;

     parameter Real G= (((substanceData.DfH - substanceData.DfG)/298.15)
       + (500000.* substanceData.E)/T0^2
       - (1000*substanceData.A_)/((-1 + exp((1000*substanceData.E_)/T0))*substanceData.E_*T0)
       - 0.001*substanceData.B*T0 - 5*10^(-7) * substanceData.C * T0^2
       - (1/3)*10^(-9)*substanceData.D*T0^3 - sqrt(0.004*T0)* substanceData.X
       + (substanceData.A_*log(1 - exp(-((1000*substanceData.E_)/T0))))/substanceData.E_^2
       - A*log(0.001*T0));

   algorithm
     //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
     //+ (substanceData.Cp+Modelica.Constants.R)*log(T/298.15);

       //Molar entropy:
       // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
       // - pressure shift: to reach the ideal gas equation at constant temperature Vm*dP = -T*dS (small amount of work)
       molarEntropyPure := G
         + A*log(t) + substanceData.B*t + substanceData.C*t^2/2 + substanceData.D*t^3/3
         - substanceData.E/(2*t^2)
         + 2*substanceData.X*t^0.5 + substanceData.A_/substanceData.E_/t/(exp(substanceData.E_/t) - 1)
         - substanceData.A_/substanceData.E_^2*log(1 - exp(-substanceData.E_/t))
       - Modelica.Constants.R*log(p/100000);

   /*    AA*Log[t] + BB*t + CC*t^2/2 + DD*t^3/3 - EE/(2*t^2) + 2*X*t^0.5 + G +
 AAA/EEE/t/(Exp[EEE/t] - 1) - AAA/EEE^2*Log[1 - Exp[-EEE/t]]

 G + AA*Log[t] + BB*t + CC*t^2/2 + DD*t^3/3 - EE/(2*t^2) + 2*X*t^0.5 +
 AAA/EEE/t/(Exp[EEE/t] - 1) - AAA/EEE^2*Log[1 - Exp[-EEE/t]]
 */

       //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
       //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
   end molarEntropyPure;

   redeclare function extends molarVolumePure
    "Molar volume of the pure substance"
   algorithm
       molarVolumePure := Modelica.Constants.R*T/p; //ideal gas
   end molarVolumePure;

   redeclare function extends molarHeatCapacityCp
    "Molar heat capacity of the substance at constant pressure"
    protected
     parameter Real T0=298.15;
     Real t=T/1000;
     parameter Real A= substanceData.Cp
       - ((10^6 * substanceData.A_* exp(1000*substanceData.E_)/T0)) / ((-1 + exp((1000*substanceData.E_)/T0))^2 * T0^2)
       - (10^6 * substanceData.E)/T0^2 - 0.001*substanceData.B*T0 - 10^(-6) * substanceData.C * T0^2
       - 10^(-9) * substanceData.D * T0^3 - sqrt(1/1000)* T0^0.5 * substanceData.X;
   algorithm
       molarHeatCapacityCp := (A + substanceData.B*t + substanceData.C*t^2 +
       substanceData.D*t^3 + substanceData.E/t^2 + substanceData.X*t^0.5 +
       substanceData.A_/t^2*exp(substanceData.E_/t)/(exp(substanceData.E_/t)-1)^2);
   end molarHeatCapacityCp;

   redeclare function extends molarMassOfBaseMolecule "Molar mass of the substance"
   algorithm
       molarMass := substanceData.MolarWeight;
   end molarMassOfBaseMolecule;

   redeclare function extends temperature "Temperature of substance from its enthalpy"
    protected
        function f_nonlinear "Solve molarEnthalpy(data,T) for T with given molar enthalpy"
          extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
          input SubstanceData data "Ideal gas data";
          input Modelica.Units.SI.SpecificEnthalpy h "Specific enthalpy";
        algorithm
          y := specificEnthalpy(data,u)
                 - h;
        end f_nonlinear;

   algorithm
     T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
          function f_nonlinear(data=substanceData, h=h), 200, 6000);

   end temperature;

    redeclare function extends solution_temperature
    "Temperature of the solution from enthalpies os substances"
    protected
      Modelica.Units.SI.MolarMass bMM[size(X,1)] = molarMassOfBaseMolecule(substanceData);
      //this is gas, so the self-clustering is not included:
      Modelica.Units.SI.MoleFraction x[size(X,1)]=(X./molarMassOfBaseMolecule(substanceData))/
        sum(X./molarMassOfBaseMolecule(substanceData)) "mole fractions of substances";
      SubstanceData solutionData= SubstanceData(
             MolarWeight = sum(x[i]*substanceData[i].MolarWeight for i in 1:size(X,1)),
             z = sum(x[i]*substanceData[i].z for i in 1:size(X,1)),
             DfG = sum(x[i]*substanceData[i].DfG for i in 1:size(X,1)),
             DfH = sum(x[i]*substanceData[i].DfH for i in 1:size(X,1)),
             gamma =  sum(x[i]*substanceData[i].gamma for i in 1:size(X,1)),
             Cp =  sum(x[i]*substanceData[i].cp_25degC for i in 1:size(X,1)),
             B =  sum(x[i]*substanceData[i].B for i in 1:size(X,1)),
             C =  sum(x[i]*substanceData[i].C for i in 1:size(X,1)),
             D =  sum(x[i]*substanceData[i].D for i in 1:size(X,1)),
             E =  sum(x[i]*substanceData[i].E for i in 1:size(X,1)),
             X =  sum(x[i]*substanceData[i].X for i in 1:size(X,1)),
             A_ = sum(x[i]*substanceData[i].A_ for i in 1:size(X,1)),
             E_ = sum(x[i]*substanceData[i].E_ for i in 1:size(X,1)));      //TODO: gamma,X,E_ are only estimations
    algorithm
      assert(abs(sum(X))<1e-5,"sum(X) must be 1");
      T := temperature(solutionData,h,p,v,I);
    end solution_temperature;

   redeclare function extends density
       "Return density of the substance in the solution"
   algorithm
          density := substanceData.MolarWeight * p / (Modelica.Constants.R * T);
   end density;

    annotation (Documentation(revisions="<html>
<p><i>2016</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end IdealGasShomate;

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

  model Total "Summation of all extensible properties per substance"
    replaceable package stateOfMatter =
        Chemical.Interfaces.StateOfMatter constrainedby StateOfMatter
    "Substance model to translate data into substance properties"
       annotation (choicesAllMatching = true);

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

   replaceable package stateOfMatter =
        Incompressible
      constrainedby StateOfMatter
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

    parameter Boolean ElectricGround = true
    "Is electric potential equal to zero?"
      annotation (Evaluate=true, choices(checkBox=true), Dialog(group="Environment relationships"));

  Modelica.Units.SI.Temperature temperature "Temperature";

  Modelica.Units.SI.Pressure pressure "Pressure";

  Modelica.Units.SI.Volume volume "Current volume of the solution";

  Modelica.Units.SI.Mass mass(stateSelect=StateSelect.prefer)
    "Current mass of the solution";

    Total total(redeclare package stateOfMatter =
          stateOfMatter, ElectricGround=ElectricGround)
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

    extends Interfaces.PartialSolution;

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
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),
            Dialog(__Dymola_compact=true));

  parameter Modelica.Units.SI.MolarFlowRate SubstanceFlow=0
    "Volumetric flow of Substance if useSubstanceFlowInput=false"
    annotation (HideResult=true, Dialog(enable=not
          useSubstanceFlowInput));

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

  partial model ConditionalKinetics
    "Input of kinetics coefficient vs. parametric kinetics coefficient"

    parameter Boolean useForwardRateInput=false
      "= true, if forward rate coefficient is provided via input"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),
        Dialog(group="Chemical kinetics", __Dymola_compact=true));

    parameter Real k_forward(unit="mol/s") "Forward rate coefficient (mole-fraction based)  if useForwardRateInput=false"
      annotation (HideResult=true, Dialog(group="Chemical kinetics", enable=not useForwardRateInput));

    Modelica.Blocks.Interfaces.RealInput kfInput(start=KC, final unit="mol2.s-1.J-1")=
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

  connector Inlet "Electro-chemical potential and molar change of the substance in the solution"

  Modelica.Units.SI.ChemicalPotential r
    "Inertial Electro-chemical potential";

  flow Modelica.Units.SI.MolarFlowRate n_flow
    "Molar change of the substance";

  input Chemical.Utilities.Units.URT uRT "u/(R*T)";
    // u .. electro-chemical potential of the substance in solution
    // R .. gas constant
    // T .. temperature


  input Modelica.Units.SI.MolarEnthalpy h
    "Enthalphy of the substance";



    annotation (Icon(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{-100,100},{-40,0},{-100,-100},{100,0},{-100,100}},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5,
            lineColor={158,66,200})}),
      Diagram(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{52,0},{-48,50},{-28,0},{-48,-50},{52,0}},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5,
            lineColor={158,66,200})}),
      Documentation(revisions="<html>
<p><i>2023</i></p>
<p>Marek Matejak </p>
</html>",   info="<html>

<p>Chemical streams:</p>
<h4>u = û + r</h4>
<h4>r = der(q)*L</h4>
<p>u .. electro-chemical potential</p>
<p>û .. steady-state electro-chemical potential</p>
<p>r .. electro-chemical inertia</p>
<p>q .. molar flow rate</p>
<p>L .. electro-chemical inductance</p>

<p>Definition of electro-chemical potential of the substance:</p>
<h4>u(x,T,v) = u&deg;(T) + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. eletric potential of the solution</p>
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
  end Inlet;

  connector Outlet "Electro-chemical potential and molar change of the substance in the solution"

  Modelica.Units.SI.ChemicalPotential r
    "Inertial Electro-chemical potential";

  flow Modelica.Units.SI.MolarFlowRate n_flow
    "Molar change of the substance";

  output Chemical.Utilities.Units.URT uRT "u/(R*T)";
    // u .. electro-chemical potential of the substance in solution
    // R .. gas constant
    // T .. temperature


  output Modelica.Units.SI.MolarEnthalpy h
    "Enthalphy of the substance";



    annotation (
      Icon(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{100,0},{-100,100},{-40,0},{-100,-100},{100,0}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={158,66,200},
            lineThickness=0.5)}),
      Diagram(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{50,0},{-50,50},{-30,0},{-50,-50},{50,0}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={158,66,200},
            lineThickness=0.5)}),
       Documentation(revisions="<html>
<p><i>2023</i></p>
<p>Marek Matejak </p>
</html>",   info="<html>

<p>Chemical streams:</p>
<h4>u = û + r</h4>
<h4>r = der(q)*L</h4>
<p>u .. electro-chemical potential</p>
<p>û .. steady-state electro-chemical potential</p>
<p>r .. inertial electro-chemical potential</p>
<p>q .. molar flow rate</p>
<p>L .. electro-chemical inductance</p>

<p>Definition of electro-chemical potential of the substance:</p>
<h4>u(x,T,v) = u&deg;(T) + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. eletric potential of the solution</p>
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
  end Outlet;

  partial model SISOFlow "Base Model with basic flow eqautions for SISO"
    import Chemical;
    import Chemical.Utilities.Types.InitializationMethods;

    parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
      annotation(Dialog(tab="Advanced"));
    parameter InitializationMethods initN_flow = Chemical.Utilities.Types.InitializationMethods.none "Initialization method for n_flow"
      annotation(Dialog(tab= "Initialization", group="Molar flow"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_0 = 0 "Initial value for n_flow"
      annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.state)));
    parameter Utilities.Units.MolarFlowAcceleration n_acceleration_0 = 0 "Initial value for der(n_flow)"
      annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.derivative)));

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the molar flow" annotation (Dialog(tab="Advanced"));

    InletProcess inlet annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
    Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(extent={{80,-20},{120,20}})));

    Modelica.Units.SI.MolarFlowRate n_flow(stateSelect=n_flowStateSelect) = inlet.n_flow
        "Molar flow through component";

    // inlet state quantities
  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Chemical.Utilities.Units.URT uRT_in(unit="1")=inlet.uRT "Electro-chemical potential of substance entering divided by (R*T)";
    Chemical.Utilities.Units.URT u0RT_in(unit="1")=inlet.u0RT "Electro-chemical potential of pure substance entering divided by (R*T)";
    Modelica.Units.SI.MolarEnthalpy h_in=inlet.h "Enthalpy of substance enetering";

    //outlet state quantities
    Chemical.Utilities.Units.URT uRT_out(unit="1") "Electro-chemical potential of substance exiting divided by (R*T)";
    Modelica.Units.SI.MolarEnthalpy h_out "Enthalpy of substance exiting";

  initial equation
    if initN_flow == InitializationMethods.state then
      n_flow = n_flow_0;
    elseif initN_flow == InitializationMethods.derivative then
      der(n_flow) = n_acceleration_0;
    elseif initN_flow == InitializationMethods.steadyState then
      der(n_flow) = 0;
    end if;
  equation

    inlet.n_flow + outlet.n_flow = 0;
    outlet.r = inlet.r - der(inlet.n_flow) * L;


    outlet.uRT = uRT_out;
    outlet.h = h_out;

    h_out = h_in;

    annotation (Documentation(info="<html>
<p>Interface class for all components with an Inlet and an Outlet and a molarflow without a mass storage between.</p>
<p>This class already implements the equations that are common for such components, namly the conservation of mass, the intertance equation. </p>
</html>"));
  end SISOFlow;

  partial model SISOFlowVertical "Base Model with basic flow eqautions for SISO"
    import Chemical;
    import Chemical.Utilities.Types.InitializationMethods;

    parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
      annotation(Dialog(tab="Advanced"));
    parameter InitializationMethods initN_flow = Chemical.Utilities.Types.InitializationMethods.none "Initialization method for n_flow"
      annotation(Dialog(tab= "Initialization", group="Molar flow"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_0 = 0 "Initial value for n_flow"
      annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.state)));
    parameter Utilities.Units.MolarFlowAcceleration n_acceleration_0 = 0 "Initial value for der(n_flow)"
      annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.derivative)));

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the molar flow" annotation (Dialog(tab="Advanced"));

    InletProcess inlet annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,100})));
    Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-100})));

    Modelica.Units.SI.MolarFlowRate n_flow(stateSelect=n_flowStateSelect) = inlet.n_flow
        "Molar flow through component";


    // inlet state quantities
  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Chemical.Utilities.Units.URT uRT_in(unit="1")=inlet.uRT "Electro-chemical potential of substance entering divided by (R*T)";
    Chemical.Utilities.Units.URT u0RT_in(unit="1")=inlet.u0RT "Electro-chemical potential of pure substance entering divided by (R*T)";
    Modelica.Units.SI.MolarEnthalpy h_in=inlet.h "Enthalpy of substance enetering";

    //outlet state quantities
    Chemical.Utilities.Units.URT uRT_out(unit="1") "Electro-chemical potential of substance exiting divided by (R*T)";
    Modelica.Units.SI.MolarEnthalpy h_out "Enthalpy of substance exiting";

  initial equation
    if initN_flow == InitializationMethods.state then
      n_flow = n_flow_0;
    elseif initN_flow == InitializationMethods.derivative then
      der(n_flow) = n_acceleration_0;
    elseif initN_flow == InitializationMethods.steadyState then
      der(n_flow) = 0;
    end if;
  equation

    inlet.n_flow + outlet.n_flow = 0;
    outlet.r = inlet.r - der(inlet.n_flow) * L;

    outlet.uRT = uRT_out;
    outlet.h = h_out;

    h_out = h_in;

    annotation (Documentation(info="<html>
<p>Interface class for all components with an Inlet and an Outlet and a molarflow without a mass storage between.</p>
<p>This class already implements the equations that are common for such components, namly the conservation of mass, the intertance equation. </p>
</html>"));
  end SISOFlowVertical;

  connector InletProcess "Inlet with formation energy of the substance"

  Modelica.Units.SI.ChemicalPotential r
    "Inertial Electro-chemical potential";

  flow Modelica.Units.SI.MolarFlowRate n_flow
    "Molar change of the substance";

  input Chemical.Utilities.Units.URT uRT "u/(R*T)";
    // u .. electro-chemical potential of the substance in solution
    // R .. gas constant
    // T .. temperature

  input Modelica.Units.SI.MolarEnthalpy h
    "Enthalphy of the substance";

  input Chemical.Utilities.Units.URT u0RT "u0/(R*T)";
    // u0 .. electro-chemical potential of the pure substance
    // R .. gas constant
    // T .. temperature

    annotation (Icon(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{-100,100},{-40,0},{-100,-100},{100,0},{-100,100}},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5,
            lineColor={158,66,200})}),
      Diagram(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{52,0},{-48,50},{-28,0},{-48,-50},{52,0}},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5,
            lineColor={158,66,200})}),
      Documentation(revisions="<html>
<p><i>2023</i></p>
<p>Marek Matejak </p>
</html>",   info="<html>

<p>Chemical streams:</p>
<h4>u = û + r</h4>
<h4>r = der(q)*L</h4>
<p>u .. electro-chemical potential</p>
<p>û .. steady-state electro-chemical potential</p>
<p>r .. electro-chemical inertia</p>
<p>q .. molar flow rate</p>
<p>L .. electro-chemical inductance</p>

<p>Definition of electro-chemical potential of the substance:</p>
<h4>u(x,T,v) = u&deg;(T) + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. eletric potential of the solution</p>
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
  end InletProcess;

  connector OutletSubstance
    "Outlet with formation energy of the substance"

  Modelica.Units.SI.ChemicalPotential r
    "Inertial Electro-chemical potential";

  flow Modelica.Units.SI.MolarFlowRate n_flow
    "Molar change of the substance";

  output Chemical.Utilities.Units.URT uRT "u/(R*T)";
    // u .. electro-chemical potential of the substance in solution
    // R .. gas constant
    // T .. temperature


  output Modelica.Units.SI.MolarEnthalpy h
    "Enthalphy of the substance";


  output Chemical.Utilities.Units.URT u0RT "u0/(R*T)";
    // u0 .. electro-chemical potential of the pure substance
    // R .. gas constant
    // T .. temperature

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{100,0},{-100,100},{-40,0},{-100,-100},{100,0}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={158,66,200},
            lineThickness=0.5)}),
      Diagram(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{50,0},{-50,50},{-30,0},{-50,-50},{50,0}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={158,66,200},
            lineThickness=0.5)}),
       Documentation(revisions="<html>
<p><i>2023</i></p>
<p>Marek Matejak </p>
</html>",   info="<html>

<p>Chemical streams:</p>
<h4>u = û + r</h4>
<h4>r = der(q)*L</h4>
<p>u .. electro-chemical potential</p>
<p>û .. steady-state electro-chemical potential</p>
<p>r .. inertial electro-chemical potential</p>
<p>q .. molar flow rate</p>
<p>L .. electro-chemical inductance</p>

<p>Definition of electro-chemical potential of the substance:</p>
<h4>u(x,T,v) = u&deg;(T) + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. eletric potential of the solution</p>
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
  end OutletSubstance;
end Interfaces;
