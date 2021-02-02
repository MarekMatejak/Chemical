within ;
package Chemical "Physical Chemistry"
  package UsersGuide "User's Guide"
    extends Modelica.Icons.Information;

  class Overview "Overview"
    extends Modelica.Icons.Information;

   annotation (Documentation(info="<html>
<p>The Chemical library can describe the following phenomena.</p>
<table cellspacing=\"0\" cellpadding=\"2\" border=\"1\"><tr>
<td><h4>Chemical Components</h4></td>
<td><h4>Description</h4></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Solution1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Solution\">Chemical solution</a></p><p>The solution is the base component of each model, because it defines the conditions of the electro-chemical processes. It integrates the total amount of substance (called amount of solution), heat, charge, entropy, volume and others from each substances to present the base properties such as temperature, pressure, electric potential and others. The usage is very simple - just connect each chemical substance with its chemical solution using their <a href=\"modelica://Chemical.Interfaces.SolutionPort\">SolutionPort</a>.</p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Substance1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Substance\">Chemical substance</a></p><p>The chemical substance integrates the amount of the chemical substance and from the properties of the connected solution it presents the electro-chemical potential of the substance using the <a href=\"modelica://Chemical.Interfaces.SubstancePort\">SubstancePort</a>.</p><p>There are two basic <a href=\"modelica://Chemical.Interfaces.StateOfMatter\">states of matter</a>: <a href=\"modelica://Chemical.Interfaces.IdealGas\">ideal gas</a> and <a href=\"modelica://Chemical.Interfaces.Incompressible\">incompressible</a> substance. However, the user can easily (re)define their own state of matter by inserting the correct expressions for the pure substance <a href=\"modelica://Chemical.Interfaces.StateOfMatter.activityCoefficient\">activity coefficient</a>, <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarVolumePure\">molar volume</a>, <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarEntropyPure\">molar entropy</a> and <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarEnthalpyElectroneutral\">molar enthalpy</a>, based on the current solution state (temperature, pressure, electric potential and ionic strength) and the <a href=\"modelica://Chemical.Interfaces.StateOfMatter.SubstanceData\">substance data</a>. The object-oriented design allows users to define the substance data record as part of the state of matter package. Users can select substance parameters according to the state of matter, redefining the getter functions of substance properties.</p><p>The examples work with ideal gases in case of all gaseous substance and incompressible state of matter in case of liquid or solid. The definition data are the molar mass of the substance, the number of charges of the substance, the molar heat capacity of the substance at a constant pressure, free formation enthalpy, free formation Gibbs energy and density (if incompressible) &mdash; all at a temperature of 25&deg;C and pressure 1 bar. Since these parameters are usually recorded in chemical tables at this standard conditions. In this manner, more than 35 real chemical <a href=\"modelica://Chemical.Examples.Substances\">substances</a> in the example package of this chemical library have already been defined. The usage of these predefined substances&rsquo; data is very simple. In the parameter dialog of the chemical substance, the correct record with this data can be selected, as shown in Figure 1.</p><p>This setting is typically the most important setting of each chemical model. All equilibrium coefficients, standard voltages, dissolution coefficients, saturated vapor pressures and so on, are automatically solved using these substance data. As a result, for example, the chemical reaction component only needs to define the stoichiometry coefficients, and the connected substances reach equilibrium at the correct equilibrium coefficient.</p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Reaction1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Reaction\">Chemical reaction</a></p><p>The chemical reaction component is very general. The dissociation constant of the equilibrium is calculated from substance properties at usual in thermodynamics, for example as definition of <a href=\"http://goldbook.iupac.org/S05915.html\">UIPAC</a>. For example if we want to define <a href=\"modelica://Chemical.Examples.SimpleReaction\">simple reaction A&lt;-&gt;B</a> with dissociation constant [B]/[A]=2 then it must be the difference between Gibbs energies of formation equal to B.DfG - A.DfG = - R * T * ln(2). Without lost of generality it is possible to select some substances as reference and give them the zero Gibbs energy of formation. The next substances created by some chemical process can be expressed from them such as example of <a href=\"modelica://Chemical.Examples.Hemoglobin.Allosteric_Hemoglobin_MWC\">alosteric hemoglobin</a> calculation. The kinetics of the chemical reaction is different as usual. However the most of processes can be recalculated with sufficient precision, for example the <a href=\"Chemical.Examples.MichaelisMenten\">Michaelic-Menton</a> can be recalculated with precision of 1.5% of maximal rate. </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Diffusion1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Diffusion\">Diffusion</a></p><p>Diffusion is a dynamic chemical process, wich is also equilibrating of electro-chemical potential of the substance. Analogically as in chemical reaction the speed of diffucion can be calculated as coefficient C multiplied by electro-chemical gratient. C can be a parammeter or input expressed from distance, substance and solution properties. </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/GasSolubility1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.GasSolubility\">Henry&apos;s law, Raoult&apos;s law or Sieverts&apos; law</a></p><p>Surprisingly, all these laws has the same basis = equilibrium of electro-chemical potential. The most of problems in data is caused by wrong selection of standard state as 1 mol/kg or 1 mol/L. Please avoid these assumptions of these totally confused states and use only mole fractions instead of each molality or molarity - the world will be much better (I promise). </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Membrane1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Membrane\">Semipermeable membrane</a></p><p>The same as before - just equilibrating the electro-chemical potentials. A result is the Donnan&apos;s equilibrium, Nernst potentials of the ions and the membrane electric potential. Transporting water through membrane is reaching the osmotic equilibrium (The real one, not the simplified one defined by osmotic pressure lineary dependent on impermeable substance concentration). </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Speciation1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Speciation\">Chemical speciation</a></p><p>The chemical speciation is for macromolecule composed with independent subunits is specific conformations. For example the hemoglobin is tetramer, which can be in two conformation: relaxed and tensed. In each of this conformation it has different afinities (different dissociation constant) for binding oxygen in each of four independent subunits. This alosteric effect can be modeled using speciation such as in <a href=\"modelica://Chemical.Examples.Hemoglobin.Allosteric_Hemoglobin2_MWC\">Allosteric_Hemoglobin2_MWC</a>. However the result should be the same as using the detailed reaction model <a href=\"modelica://Chemical.Examples.Hemoglobin.Allosteric_Hemoglobin_MWC\">Allosteric_Hemoglobin_MWC</a>.</p></td>
</tr>
</table>
</html>"));
  end Overview;

  class Connectors "Connectors"
    extends Modelica.Icons.Information;

   annotation (Documentation(info="<html>
<p>The Chemical defines the two important <b>elementary connectors</b> for substance and for solution:</p>
<table cellspacing=\"0\" cellpadding=\"1\" border=\"1\"><tr>
<td valign=\"top\"></td>
<td valign=\"top\"><h4>potential</h4><p>variables</p></td>
<td valign=\"top\"><h4>flow</h4><p>variables</p></td>
<td valign=\"top\"><h4>stream</h4><p>variables</p></td>
<td valign=\"top\"><h4>connector definition</h4></td>
<td valign=\"top\"><h4>icons</h4></td>
</tr>
<tr>
<td valign=\"middle\"><h4>substance</h4></td>
<td valign=\"middle\"><p>u .. electro-chemical potential of the chemical substance</p></td>
<td valign=\"middle\"><p>q .. molar flow of the chemical substance</p></td>
<td valign=\"middle\"></td>
<td valign=\"middle\"><p><br><a href=\"Chemical.Interfaces.SubstancePort\">Chemical.Interfaces.SubstancePort</a> </p></td>
<td valign=\"middle\"><p><img src=\"modelica://Chemical/Resources/Images/UsersGuide/ChemicalPorts.png\"/></p></td>
</tr>
<tr>
<td valign=\"middle\"><h4>solution</h4></td>
<td valign=\"middle\"><p>p .. pressure of the solution</p><p>T .. temperature of the solution</p><p>v .. electric potential of the solution</p><p><br>n .. amount of all substances in the solution</p><p>m .. mass of the solution</p><p>V .. volume of the solution</p><p>G .. free Gibbs energy of the solution</p><p>Q .. electric charge of the solution</p><p>I .. ionic strength of the solution</p></td>
<td valign=\"middle\"><p>dV .. change of the volume of the solution</p><p>dH .. enthalpy change of the solution</p><p>i .. electric charge change of the solution</p><p><br><i>nj ..  amount of the substance</i></p><p><i>mj .. mass of the substance</i></p><p><i>Vj .. volume of the substance</i></p><p><i>Gj .. free Gibbs energy of the substance</i></p><p><i>Qj .. electric charge of the substance</i></p><p><i>Ij .. ionic strength of the substance</i></p></td>
<td valign=\"middle\"></td>
<td valign=\"middle\"><p><br><a href=\"Chemical.Interfaces.SolutionPort\">Chemical.Interfaces.SolutionPort</a></p></td>
<td valign=\"middle\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/SolutionPort.png\"/></p></td>
</tr>
<tr>
<td valign=\"middle\"><h4>substanceMass</h4></td>
<td valign=\"middle\"><p>x_mass .. mass fraction of the chemical substance in solution</p></td>
<td valign=\"middle\"><p>m_flow .. mass flow of the chemical substance</p></td>
<td valign=\"middle\"></td>
<td valign=\"middle\"><p><br><a href=\"Chemical.Interfaces.SubstanceMassPort\">Chemical.Interfaces.SubstanceMassPort</a> </p></td>
<td valign=\"middle\"><p><img src=\"modelica://Chemical/Resources/Images/UsersGuide/ChemicalMassPorts.png\"/></p></td>
</tr>
<tr>
<td valign=\"middle\"><h4>substanceMolarity</h4></td>
<td valign=\"middle\"><p>c .. molar concentration per liter of the chemical substance in solution</p></td>
<td valign=\"middle\"><p>q .. molar flow of the chemical substance</p></td>
<td valign=\"middle\"></td>
<td valign=\"middle\"><p><br><a href=\"Chemical.Interfaces.SubstanceMolarityPort\">Chemical.Interfaces.SubstanceMolarityPort</a> </p></td>
<td valign=\"middle\"><p><img src=\"modelica://Chemical/Resources/Images/UsersGuide/ChemicalMolarityPorts.png\"/></p></td>
</tr>
</table>
</html>"));
  end Connectors;

  package ReleaseNotes "Release notes"
    extends Modelica.Icons.ReleaseNotes;

  class Version_1_0 "Version 1.0.0 (Apr. 28, 2015)"
    extends Modelica.Icons.ReleaseNotes;

  annotation (Documentation(info="<html>
<ul>
<li>Separation the Chemical from Physiolibrary to https://github.com/MarekMatejak/Chemical from https://github.com/MarekMatejak/Physiolibrary branche PhysicalChemistry </li>
<li><font style=\"color: #333333; \">Components for solution, substance, chemical reaction, diffusion, gas dissolution, semipermeable membranes, chemical speciation of macromolecules, ..</font></li>
<li><font style=\"color: #333333; \">The library uses the Modelica Standard Libary (MSL) version 3.2.</font></li>
</ul>
</html>"));
  end Version_1_0;

  class Version_1_1 "Version 1.1.0 (Sep. 15, 2015)"
    extends Modelica.Icons.ReleaseNotes;

  annotation (Documentation(info="<html>
<ul>
<li>New state of matter - ideal gas</li>
<li>Solution with mechanical, thermal and electrical port for any state of matter</li>
<li>Sensor of partial pressure</li>
<li>Sensor of dissociation constant of chemical reaction for hypothetical pure substances&rsquo; scheme </li>
<li>New Examples </li>
<li>New icon for electron transfer</li>
<li>New icon for chemical buffer</li>
<li>New chemical kinetics with speed turnover</li>
</ul>
</html>"));
  end Version_1_1;

  class Version_1_2 "Version 1.2.0 (Oct. 15, 2018)"
    extends Modelica.Icons.ReleaseNotes;

  annotation (Documentation(info="<html>
<ul>
<li>Substance, which are making clusters (E.g. liquid water molecules with hydrogen bonds</li>
<li>Support of fluid connectors</li>
<li>Mass fraction connector</li>
<li>Molarity concentration conector</li>
</ul>
</html>"));
  end Version_1_2;

  class Version_1_3 "Version 1.3.0 (Nov. 19, 2020)"
    extends Modelica.Icons.ReleaseNotes;

  annotation (Documentation(info="<html>
<ul>
<li>Elastic compartment solution (prepared for new Physiolibrary)</li>
<li>Enthalpy streams in substance connector</li>
<li>Methanogenesis examples</li>
</ul>
</html>"));
  end Version_1_3;

  class Version_1_4 "Version 1.4.0 (Jan. 27, 2021)"
    extends Modelica.Icons.ReleaseNotes;

  annotation (Documentation(info="<html>
<ul>
<li>1241 gaseous substances referenced to MSL data (Chemical.Substances.IdelGasesMSL)</li>
<li>fix of enthalpy streams</li>
<li>fix vaporization with water clustering</li>
<li>specific properties calculation for self-clustering substances</li>
<li>fix usage of molar mass (e.g. as mass of base molecule for self clustering substance)</li>
<li>changed of some function interfaces for StateOfMatter</li>
</ul>
</html>"));
  end Version_1_4;
   annotation (Documentation(info="<html>
<p>This section summarizes the changes that have been performed on the Chemical. </p>
</html>"));

  end ReleaseNotes;

  class Contact "Contact"
    extends Modelica.Icons.Contact;

   annotation (Documentation(info="<html>
<p>Marek Matej&aacute;k</p>
<p>email: marek@matfy.cz</p>
<p>skype: marek.matejak</p>
<p>tel: +420 776 301 395</p>
</html>"));

  end Contact;

    class License "BSD 3-Clause License"
       extends Modelica.Icons.Information;
      annotation (Documentation(info="<html>
<p>All files in this directory (Physiolibrary) and in all subdirectories, especially all files that build package &quot;Physiolibrary&quot; are licensed by <u><b>Marek Matejak</b></u> under the <a href=\"http://opensource.org/licenses/BSD-3-Clause\">BSD 3-Clause License</a> (with exception of files &quot;Resources/*&quot;). </p>
<h4>Licensor:</h4>
<p>Marek Matej&aacute;k,</p>
<p>Hviezdoslavova 632/41,</p>
<p>916 01 Star&aacute; Tur&aacute;, </p>
<p>Slovak Republic,</p>
<p>email: marek@matfyz.cz</p>
<h4><span style=\"color:#008000\">Organization: </span></h4>
<p>Institute of Pathological Physiology, First Faculty of Medicine, Charles University in Prague,</p>
<p>U Nemocnice 5, 128 53 Prague 2, Czech Republic</p>
<br><h4>Copyright notices of the files:</h4>
<p>Copyright (c) 2008-2015, Marek Matej&aacute;k, Charles University in Prague</p>
<p><br>All rights reserved. </p>
<p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: </p>
<p>1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. </p>
<p>2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. </p>
<p>3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. </p>
<p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS &quot;AS IS&quot; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
</html>"));
    end License;

  class NewRelease "Publishing new release"
    extends Modelica.Icons.Information;

   annotation (Documentation(info="<html>
<p><br>New release must be numbered by Semantic Versioning 2.0.0, see <a href=\"http://semver.org/\">semver.org</a>. </p>
<p><br>If minor version, then the conversion script must be written and connected with package Chemical using &quot;annotation(conversion(from(version=..)))&quot;! </p>
<p><br>To clean the code from dummy annotations try to use script <a href=\"https://github.com/dietmarw/trimtrailingwhitespaces\">ttws</a>. </p>
<p>To check english spelling try to use <a href=\"https://github.com/vlajos/misspell_fixer\">missspell_fixer</a>.</p>
<p><br>Update version number to &quot;X.Y.Z&quot;: </p>
<ul>
<li>At package Chemical annotation: (version=&quot;X.Y.Z&quot;) together with &quot;versionBuild&quot;, &quot;versionDate&quot; and &quot;dateModified&quot; attribute </li>
<li>At file &quot;./Chemical/libraryinfo.mos&quot; </li>
</ul>
<p><br>Update release notes: </p>
<ul>
<li>At UsersGuide.ReleaseNotes</li>
<li>At file &quot;./README.md&quot;, together with update of &quot;Current release&quot; section.</li>
</ul>
<p><br>Publish release in GitHub: </p>
<ul>
<li>Prepare release in &quot;master&quot; branch</li>
<li>Install, Check, Test, Test, Test.. </li>
<li>Draft a new <a href=\"https://github.com/impact/impact/blob/master/resources/docs/modelica2015/paper/impact.md#impact-on-library-developers\">release from &quot;master&quot;</a> branch with number &quot;vX.Y.Z&quot; and with release notes. </li>
</ul>
</html>"));
  end NewRelease;
  annotation (DocumentationClass=true, Documentation(info="<html>
<p>Package <b>Chemical </b>is a modelica package for <b>Electro-Chemical processes </b>that is developed from <b>Physiolibrary</b> modelica implementation, see <a href=\"http://patf-biokyb.lf1.cuni.cz/wiki/hummod/hummod\">http://www.physiolibrary.org</a>. It provides connectors and model components fitted for electro-chemical models. </p>
</html>"));
  end UsersGuide;

 extends Modelica.Icons.Package;

  package Components "Chemical Components"
    model Solution "Chemical solution as homogenous mixture of the substances"
      extends Icons.Solution;

      extends Interfaces.PartialSolutionWithHeatPort(pressure(start=BasePressure));

    parameter Modelica.Units.SI.Pressure BasePressure=system.p_ambient
      "Pressure at zero mechanic force (or if not useMechanicPorts)"
      annotation (HideResult=true);

      parameter Boolean useMechanicPorts = false "Are mechanic ports pressent?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Area SurfaceArea=0.01
      "Area for surfacePort to connect MultiBody components"
      annotation (HideResult=true, Dialog(enable=useMechanicPorts));

      parameter Boolean isPistonPositionAbsolute=false
      "Relavite position has zero at initial state without force"
        annotation (HideResult=true,choices(checkBox=true), Dialog(enable=useMechanicPorts));

      Modelica.Mechanics.Translational.Interfaces.Flange_a surfaceFlange(f=f,s=top_s) if useMechanicPorts
      "The pressure of solution generate force on prescribed surface."
        annotation (Placement(transformation(extent={{-10,70},{10,90}}),
            iconTransformation(extent={{-2,98},{2,102}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_b bottom(f=-f,s=top_s - ds) if useMechanicPorts
      "Fix of the cilinder on bottom."   annotation (Placement(transformation(
              extent={{-10,-90},{10,-70}}), iconTransformation(extent={{-2,-104},{2,
                -100}})));

      Interfaces.SolutionPort solution "Solution nonflows and flows"
                                      annotation (Placement(
            transformation(extent={{50,-90},{70,-70}}),  iconTransformation(extent={{58,-100},
              {62,-96}})));
  protected
    parameter Modelica.Units.SI.Position positionShift(fixed=false)
      "=0 absolute, otherwise negative";
    Modelica.Units.SI.Position top_s;
    Modelica.Units.SI.Position ds;
    Modelica.Units.SI.Force f;

    initial equation
      positionShift= if
                       (isPistonPositionAbsolute) then 0 else volume/SurfaceArea;

    equation

      //hydraulic
      ds = volume/SurfaceArea - positionShift;
      pressure = BasePressure - f/SurfaceArea;

      if not useMechanicPorts then
        f=0;
        top_s=ds; //equivalent for bottom_s==0
      end if;

      connect(solution, total.solution) annotation (Line(points={{60,-80},{60,-94},
              {84,-94},{84,-86}}, color={127,127,0}));
                                                                                                        annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false, initialScale=1, extent={{-100,-100},{
              100,100}}),
            graphics={Text(
              extent={{-90,-88},{78,-96}},
              lineColor={128,0,255},
              textString="%name",
              horizontalAlignment=TextAlignment.Left)}),
        Documentation(revisions="<html>
<p>2015-2020 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<h4>amountOfSolution = &sum; amountOfSubstances</h4>
<h4>mass = &sum; massOfSubstances</h4>
<h4>volume = &sum; volumeOfSubstances</h4>
<h4>freeGibbsEnergy = &sum; freeGibbsEnergiesOfSubstances</h4>
<p>To calculate the sum of extensive substance's properties is misused the Modelica \"flow\" prefix even there are not real physical flows. </p>
</html>"));
    end Solution;

    model Substance "Substance in solution"
      extends Icons.Substance;

      Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
        "Molar concentration of particles";

      extends Interfaces.PartialSubstanceInSolutionWithAdditionalPorts;



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
      parameter Modelica.Units.SI.Mass m_start=if use_mass_start then mass_start else
        amountOfSubstance_start*molarMassOfBaseMolecule;

      parameter Modelica.Units.SI.MolarMass molarMassOfBaseMolecule = stateOfMatter.molarMassOfBaseMolecule(substanceData);

      Modelica.Units.SI.AmountOfSubstance amountOfBaseMolecules(start=
           m_start/molarMassOfBaseMolecule)
        "Amount of base molecules inside all clusters in compartment";

      Modelica.Units.SI.AmountOfSubstance amountOfFreeMolecule(start=
           m_start*stateOfMatter.specificAmountOfFreeBaseMolecule(
                                       substanceData,
                                       T=system.T_ambient,
                                       p=system.p_ambient))
        "Amount of free molecules not included inside any clusters in compartment";

      Modelica.Units.SI.AmountOfSubstance amountOfParticles(start=
           m_start*stateOfMatter.specificAmountOfParticles(
                                       substanceData,
                                       T=system.T_ambient,
                                       p=system.p_ambient))
        "Amount of particles/clusters in compartment";

      Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
          Modelica.Constants.R*solution.T))
        "Dissociation constant of hydrogen bond between base molecules";

      Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
          stateOfMatter.selfClusteringBondEnthalpy(substanceData)
        - solution.T * stateOfMatter.selfClusteringBondEntropy(substanceData)
        "Gibbs energy of hydrogen bond between H2O molecules";

      Modelica.Units.SI.AmountOfSubstance amountOfBonds
        "Amount of hydrogen bonds between molecules in compartment";

      Real logn(stateSelect=StateSelect.prefer, start=log(m_start/molarMassOfBaseMolecule))
      "Natural logarithm of the amount of base molecules in solution";


      parameter Boolean EnthalpyNotUsed=false annotation (
        Evaluate=true,
        HideResult=true,
        choices(checkBox=true),
        Dialog(tab="Advanced", group="Performance"));


    initial equation

      amountOfBaseMolecules = m_start/molarMassOfBaseMolecule;
    equation

      if stateOfMatter.selfClustering(substanceData) then

        //Liquid cluster theory - equilibrium:
        //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules
        //amountOfParticles/solution.n = x/(1-K*x);                //sum(x[i])
        //amountOfBaseMolecules/solution.n = x/((1-K*x)^2);            //sum(i*x[i])
        //amountOfHydrogenBonds/solution.n = x*x*K/((1-K*x)^2);   //sum((i-1)*x[i])

        amountOfParticles*(1 - SelfClustering_K*x) = amountOfFreeMolecule;

        //Calculation of "abs(amountOfBaseMolecules*(1 - SelfClustering_K*x)) = amountOfParticles":
        x = ((2*SelfClustering_K+solution.n/amountOfBaseMolecules) - sqrt((4*SelfClustering_K*solution.n/amountOfBaseMolecules)+(solution.n/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));


        amountOfBonds = amountOfBaseMolecules*x*SelfClustering_K;

        //TODO: may be the volume of the same number of free water molecules is different as volume of the same number of water molecules in cluster ..
        //TODO: more precise calculation of other properties


       //der(enthalpy) = solution.dH + q*actualStream(port_a.h_outflow);
       //enthalpy = molarEnthalpy*amountOfBaseMolecules + amountOfAdditionalBonds*bondEnthalpy;
        solution.dH =if (EnthalpyNotUsed) then 0 else der(molarEnthalpy)*
          amountOfBaseMolecules + q*molarEnthalpy - q*actualStream(port_a.h_outflow) + (
          if (calculateClusteringHeat) then stateOfMatter.selfClusteringBondEnthalpy(
          substanceData)*der(amountOfBonds) else 0)
                        "heat transfer from other substances in solution [J/s]";

        solution.Gj =amountOfBaseMolecules*port_a.u + amountOfBonds*SelfClustering_dG
                        "Gibbs energy of the substance";

      else

        amountOfParticles = amountOfFreeMolecule;
        amountOfBaseMolecules = amountOfFreeMolecule;
        amountOfBonds = 0;

        //der(enthalpy) = solution.dH + q*actualStream(port_a.h_outflow);
        //enthalpy = molarEnthalpy*amountOfBaseMolecules;
        solution.dH =
          if (EnthalpyNotUsed) then  0
          else    der(molarEnthalpy)*amountOfBaseMolecules + q*molarEnthalpy
                  -q*actualStream(port_a.h_outflow)
                  "heat transfer from other substances in solution [J/s]";

        solution.Gj = amountOfBaseMolecules*port_a.u "Gibbs energy of the substance [J]";

      end if;



      //The main accumulation equation is "der(amountOfBaseMolecules)=q"
      // However, the numerical solvers can handle it in form of log(n) much better. :-)
      der(logn) = (q/amountOfBaseMolecules) "accumulation of amountOfBaseMolecules=exp(logn) [mol]";
      amountOfBaseMolecules = exp(logn);

      x = amountOfFreeMolecule/solution.n "mole fraction [mol/mol]";

      c = amountOfParticles/solution.V "concentration [mol/m3]";

      //solution flows
      solution.i = Modelica.Constants.F*z*q +
          Modelica.Constants.F*der(z)*amountOfBaseMolecules "change of sunstance charge [A]";
      solution.dV = molarVolume*q + der(molarVolume)*amountOfBaseMolecules "change of substance volume [m3/s]";

      //extensive properties
      solution.nj = amountOfParticles;
      solution.mj = amountOfBaseMolecules*molarMassOfBaseMolecule;
      solution.Vj = amountOfBaseMolecules*molarVolume;
      solution.Qj = Modelica.Constants.F*amountOfBaseMolecules*z;
      solution.Ij = (1/2)*(amountOfBaseMolecules*z^2);

         annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={Text(
              extent={{-84,22},{92,64}},
              lineColor={128,0,255},
              textString="%name")}), Documentation(revisions="<html>
<p>2009-2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<h4>n = x &middot; n(solution) = &int; MolarFlow</h4>
<p>where n is amount of the substance and x is mole fraction.</p>
<p>The main class from &ldquo;Chemical&rdquo; package is called &quot;Substance&quot;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &quot;n&quot; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &quot;n(solution)&quot; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solution. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
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
    extends Modelica.Icons.Package;

    model Reaction "Chemical Reaction"
      extends Interfaces.ConditionalKinetics;

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

      parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient"
        annotation(Dialog(group="Chemical kinetics"));

      Modelica.Units.SI.MolarFlowRate rr(start=0) "Reaction molar flow rate";

      Interfaces.SubstancePorts_b substrates[nS] annotation (Placement(
          transformation(extent={{-10,-40},{10,40}},
          rotation=180,
          origin={-100,0}),                             iconTransformation(
            extent={{-10,-40},{10,40}},
          rotation=180,
          origin={-100,0})));

      Interfaces.SubstancePorts_b products[nP] annotation (Placement(
          transformation(extent={{-10,-40},{10,40}},
          rotation=180,
          origin={100,0}),                            iconTransformation(extent={{-10,-40},
              {10,40}},
          rotation=180,
          origin={100,0})));

      Modelica.Units.SI.MolarEnthalpy h_mix;

      parameter Boolean EnthalpyNotUsed=false annotation (
        Evaluate=true,
        HideResult=true,
        choices(checkBox=true),
        Dialog(tab="Advanced", group="Performance"));
  protected
      Modelica.Units.SI.ChemicalPotential du;
    equation
      //the main equation
      du = ((p * products.u) - (s * substrates.u));
      rr = - kC * du * exp(-kE*abs(du));

      //reaction molar rates
      rr*s = substrates.q;
      rr*p = -products.q;


      // Implicit definition of the inStream()operator applied to inside connector i

      substrates.h_outflow = h_mix*ones(nS);
      products.h_outflow = h_mix*ones(nP);

      if
        (rr>0 and not EnthalpyNotUsed) then
        h_mix*(products.q*ones(nP)) + substrates.q*inStream(substrates.h_outflow) = 0;
      elseif
            (rr<0 and not EnthalpyNotUsed) then
        h_mix*(substrates.q*ones(nS)) + products.q*inStream(products.h_outflow) = 0;
      else
        h_mix=0;
      end if;

      // 0 = substrates.q * actualStream(substrates.h_outflow) + products.q * actualStream(products.h_outflow);
    /*  0 = sum(substrates[j].q*(if
                             (substrates[j].q > 0) then h_mix else inStream(substrates[j].h_outflow)) for j in 1:nS)
     +sum(products[k].q * (if
                             (products[k].q > 0)   then h_mix else inStream(products[k].h_outflow)) for k in 1:nP);
*/
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
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
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

    model ElectronTransfer
    "Electron transfer from the solution to electric circuit"
      extends Icons.ElectronTransfer;
      extends Interfaces.PartialSubstanceInSolution(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible,
        final substanceData = Chemical.Interfaces.Incompressible.SubstanceData(
        MolarWeight=5.4857990946e-7,
        z=-1,
        DfH=0,
        DfG=0,
        Cp=0,
        density=1e20));

      Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
          Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
              extent={{-110,-10},{-90,10}})));

    equation

      //electric
      pin.v = electricPotential;
      pin.i + z*Modelica.Constants.F*port_a.q + solution.i = 0;

      //pure substance
      x = 1;

      //solution changes
      solution.dH = 0;
      solution.dV = 0;

      //extensive properties of the solution
      solution.nj=0;
      solution.mj=0;
      solution.Vj=0;
      solution.Gj=0;
      solution.Qj=0;
      solution.Ij=0;

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Text(
              extent={{-146,-44},{154,-84}},
              textString="%name",
              lineColor={128,0,255})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ElectronTransfer;

    model Diffusion "Solute diffusion"
      extends Icons.Diffusion;
      extends Interfaces.OnePort;
      extends Interfaces.ConditionalKinetics;

      parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  protected
    Modelica.Units.SI.ChemicalPotential du;
    equation
      //the main equation
      du = (port_b.u - port_a.u);
      port_b.q = kC * du * exp(-kE*abs(du));

       annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
    end Diffusion;

    model GasSolubility "Henry's law of gas solubility in liquid."

      extends Icons.GasSolubility;

      Interfaces.SubstancePort_b gas_port "Gaseous solution"
        annotation (Placement(transformation(extent={{-10,90},{10,110}})));

      Interfaces.SubstancePort_b liquid_port "Dissolved in liquid solution" annotation (Placement(
            transformation(extent={{-10,-110},{10,-90}}), iconTransformation(extent={{-10,-110},{10,-90}})));

      extends Interfaces.ConditionalKinetics;

      parameter Real kE(unit="mol/J") = 0 "Kinetic turnover coefficient";

      parameter Boolean EnthalpyNotUsed=false annotation (
        Evaluate=true,
        HideResult=true,
        choices(checkBox=true),
        Dialog(tab="Advanced", group="Performance"));
  protected
      Modelica.Units.SI.ChemicalPotential du;
    equation
      gas_port.q + liquid_port.q = 0;

      du = (liquid_port.u - gas_port.u);

      liquid_port.q = kC*du*exp(-kE*abs(du));


      gas_port.h_outflow = if (EnthalpyNotUsed) then 0 else inStream(liquid_port.h_outflow);
      liquid_port.h_outflow = if (EnthalpyNotUsed) then 0 else inStream(gas_port.h_outflow);

      annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Gaseuous substance dissolition in liquid (Henry&apos;s law, Raoult&apos;s law, Nernst dissolution in one). </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K<sub>H</sub> =x<sub>L</sub> / x<sub>g</sub>&nbsp;</p></td>
<td><p>Henry&apos;s coefficient, Raoult&apos;s coefficient</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>G = &Delta;<sub>f</sub>G<sub>L </sub>- &Delta;<sub>f</sub>G<sub>g </sub>= &Delta;<sub>sol</sub>H - T&middot;&Delta;<sub>sol</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K<sub>H</sub>&middot; (f<sub>L</sub> / f<sub>g</sub>)) </p></td>
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
<h4><span style=\"color:#008000\">Notations</span></h4>
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
</html>"));
    end GasSolubility;

    model Membrane
    "Passive transport of the substance through semipermeable membrane"
      extends Icons.Membrane;
      extends Interfaces.OnePort;
      extends Interfaces.ConditionalKinetics;

      parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  protected
    Modelica.Units.SI.ChemicalPotential du;
    equation
      //the main equation
      du = (port_a.u - port_b.u);
      port_a.q = kC * du * exp(-kE*abs(du));

      annotation ( Documentation(info="<html>
<p><u><b><font style=\"color: #008000; \">Filtration throught semipermeable membrane.</font></b></u></p>
<p>The penetrating particles are driven by electric and chemical gradient to reach Donnan&apos;s equilibrium.</p>
<p>If zero-flow Donnan&apos;s equilibrium is reached. </p>
</html>", revisions="<html>
<p><i>2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}), graphics={
            Text(
              extent={{-97,-12},{97,12}},
              textString="%name",
              lineColor={128,0,255},
            origin={69,2},
            rotation=90)}));
    end Membrane;

    model SubstancePump "Prescribed sunstance molar flow"
      extends Interfaces.OnePort;
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
              lineColor={128,0,255},
              origin={-10,-76},
              rotation=360,
              textString="%name")}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstancePump;

    model Speciation
    "Quaternary macromolecule form defined by all its subunits"
      extends Icons.Speciation;

      replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby
      Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
         annotation (choicesAllMatching = true);

      parameter Integer NumberOfSubunits=1
      "Number of independent subunits occurring in macromolecule";

      Interfaces.SolutionPort solution                                                              annotation (Placement(transformation(extent={{-70,
                -110},{-50,-90}}),
            iconTransformation(extent={{-70,-110},{-50,-90}})));

      Modelica.Units.SI.AmountOfSubstance nm
        "Amount of the macromolecule (all form in the conformation)";
      Modelica.Units.SI.MoleFraction xm
        "Mole fraction of the macromolecule (all form of in the conformation)";

  public
      Interfaces.SolutionPort subunitSolution "The port to connect all subunits"
        annotation (Placement(transformation(extent={{-70,92},{-50,112}}),
            iconTransformation(extent={{30,50},{50,70}})));
    Interfaces.SubstancePort_a port_a annotation (Placement(transformation(
            extent={{90,-110},{110,-90}}), iconTransformation(extent={{90,-110},
              {110,-90}})));
    Interfaces.SubstancePorts_b subunits[NumberOfSubunits]
      "Subunits of macromolecule" annotation (Placement(transformation(extent={
              {-56,-14},{-36,66}}), iconTransformation(
          extent={{-10,-40},{10,40}},
          rotation=90,
          origin={-30,102})));

      parameter Boolean EnthalpyNotUsed=false annotation (
        Evaluate=true,
        HideResult=true,
        choices(checkBox=true),
        Dialog(tab="Advanced", group="Performance"));

  protected
      Modelica.Units.SI.MolarEnthalpy h_mix;
    equation
      //amount of macromolecule (all forms in conformation)
      nm*NumberOfSubunits + subunitSolution.nj = 0;

      //change of macromolecule = change of its subunits
      subunits.q = -port_a.q * ones(NumberOfSubunits);

      //mole fraction of all forms in conformation
      xm = nm/solution.n;

      //electrochemical potential of the specific form
      port_a.u = Modelica.Constants.R*solution.T*log(xm) +
            sum(subunits.u - Modelica.Constants.R*solution.T*log(xm)
             * ones(NumberOfSubunits));

      port_a.h_outflow = h_mix;
      subunits.h_outflow = (h_mix/NumberOfSubunits)*ones(NumberOfSubunits);

       if
         (port_a.q < 0 and not EnthalpyNotUsed) then
         h_mix = inStream(subunits.h_outflow) * ones(NumberOfSubunits);
       elseif
             (port_a.q > 0 and not EnthalpyNotUsed) then
         h_mix = inStream(port_a.h_outflow);
       else
         h_mix = 0;
       end if;



      //properties from subunits
      subunitSolution.dH + solution.dH = 0;
      subunitSolution.i + solution.i = 0;
      subunitSolution.Qj + solution.Qj = 0;
      subunitSolution.Ij + solution.Ij = 0;

      //properties of macromolecule as a whole
      subunitSolution.nj + solution.nj*NumberOfSubunits = 0; //only amount of substance is necessery to express between sites' solution and real solution
      subunitSolution.mj + solution.mj = 0;
      subunitSolution.Vj + solution.Vj = 0;
      subunitSolution.Gj + solution.Gj = 0;
      subunitSolution.dV + solution.dV = 0;

      //shift global solution status to subunits
      subunitSolution.T = solution.T;
      subunitSolution.v = solution.v;
      subunitSolution.p = solution.p;
      subunitSolution.n = solution.n;
      subunitSolution.m = solution.m;
      subunitSolution.V = solution.V;
      subunitSolution.G = solution.G;
      subunitSolution.Q = solution.Q;
      subunitSolution.I = solution.I;


      annotation (defaultComponentName="macromolecule",
        Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic</p>
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
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
            graphics={                                                        Text(
              extent={{-22,-106},{220,-140}},
              lineColor={128,0,255},
              textString="%name")}));
    end Speciation;

    model Stream "Flow of whole solution"
      extends Interfaces.ConditionalSolutionFlow;

      replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby
      Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
         annotation (choicesAllMatching = true);

      parameter stateOfMatter.SubstanceData substanceData
      "Definition of the substance"
         annotation (choicesAllMatching = true);

    Interfaces.SubstancePort_b port_b annotation (Placement(transformation(
            extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},
              {-90,10}})));
      Sensors.MoleFractionSensor moleFractionSensor(
         redeclare package stateOfMatter = stateOfMatter,
         substanceData=substanceData)
        annotation (Placement(transformation(extent={{56,-10},{76,10}})));
      Sensors.MoleFractionSensor moleFractionSensor1(
         redeclare package stateOfMatter = stateOfMatter,
         substanceData=substanceData)
        annotation (Placement(transformation(extent={{-56,-10},{-76,10}})));
      SubstancePump substancePump(useSubstanceFlowInput=true,EnthalpyNotUsed=EnthalpyNotUsed)
        annotation (Placement(transformation(extent={{-14,-74},{6,-54}})));
      Modelica.Blocks.Logical.Switch switch1 annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={0,-38})));
      Modelica.Blocks.Logical.GreaterThreshold greaterThreshold annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-4})));
      Modelica.Blocks.Math.Product product
        annotation (Placement(transformation(extent={{-40,-36},{-20,-16}})));
      Modelica.Blocks.Math.Product product1 annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={30,-26})));
    Interfaces.SubstancePort_b port_a
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      Interfaces.SolutionPort solution
        annotation (Placement(transformation(extent={{-70,-110},{-50,-90}})));

     parameter Boolean EnthalpyNotUsed=false annotation (
        Evaluate=true,
        HideResult=true,
        choices(checkBox=true),
        Dialog(tab="Advanced", group="Performance"));

    equation
      product.u1=q;
      product1.u2=q;
      greaterThreshold.u=q;

      connect(port_b, moleFractionSensor1.port_a) annotation (Line(
          points={{-100,0},{-76,0}},
          color={158,66,200}));
      connect(moleFractionSensor.port_a, port_a) annotation (Line(
          points={{76,0},{100,0}},
          color={158,66,200}));
      connect(moleFractionSensor1.solution, solution) annotation (Line(
          points={{-60,-10},{-60,-100}},
          color={0,128,255}));
      connect(solution, moleFractionSensor.solution) annotation (Line(
          points={{-60,-100},{60,-100},{60,-10}},
          color={0,128,255}));
      connect(substancePump.substanceFlow, switch1.y) annotation (Line(
          points={{0,-60},{0,-49},{-2.22045e-015,-49}},
          color={0,0,127}));
      connect(switch1.u2, greaterThreshold.y) annotation (Line(
          points={{2.22045e-015,-26},{0,-26},{0,-15}},
          color={255,0,255}));
      connect(product1.u1, moleFractionSensor.moleFraction) annotation (Line(
          points={{42,-32},{50,-32},{50,0},{56,0}},
          color={0,0,127}));
      connect(product.u2, moleFractionSensor1.moleFraction) annotation (Line(
          points={{-42,-32},{-50,-32},{-50,0},{-56,0}},
          color={0,0,127}));
      connect(port_b, substancePump.port_a) annotation (Line(
          points={{-100,0},{-86,0},{-86,-64},{-14,-64}},
          color={158,66,200}));
      connect(substancePump.port_b, port_a) annotation (Line(
          points={{6,-64},{84,-64},{84,0},{100,0}},
          color={158,66,200}));
      connect(product.y, switch1.u1) annotation (Line(
          points={{-19,-26},{-8,-26}},
          color={0,0,127}));
      connect(product1.y, switch1.u3) annotation (Line(
          points={{19,-26},{8,-26}},
          color={0,0,127}));
     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                            graphics={
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
              lineColor={128,0,255},
              origin={2,-74},
              rotation=180)}),
        Documentation(revisions="<html>
<p><i>2009-2018 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<h4><font color=\"#008000\">Bidirectional mass flow by concentration</font></h4>
<p>Possible field values: </p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0.1\"><tr>
<td></td>
<td><h4>forward flow</h4></td>
<td><h4>backward flow</h4></td>
</tr>
<tr>
<td><h4>solutionFlow</h4></td>
<td><p align=\"center\">&gt;=0</p></td>
<td><p align=\"center\">&lt;=0</p></td>
</tr>
<tr>
<td><h4>q_in.q</h4></td>
<td><p align=\"center\">=solutionFlow*q_in.conc</p></td>
<td><p align=\"center\">=-q_out.q</p></td>
</tr>
<tr>
<td><h4>q_out.q</h4></td>
<td><p align=\"center\">=-q_in.q</p></td>
<td><p align=\"center\">=solutionFlow*q_out.conc</p></td>
</tr>
</table>
<br/>
</html>"));
    end Stream;
  end Components;

  package Sensors "Chemical sensors"
    extends Modelica.Icons.SensorsPackage;

    model MolarFlowSensor "Measure of molar flow"

      extends Modelica.Icons.RoundSensor;
      extends Interfaces.OnePort;

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
</html>"),        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={
            Line(
              points={{70,-10},{90,-10}},
              color={127,0,127}),
            Line(
              points={{70,10},{90,10}},
              color={127,0,127}),
            Line(
              points={{-90,10},{-70,10}},
              color={127,0,127}),
            Line(
              points={{-90,-10},{-70,-10}},
              color={127,0,127}),
            Text(
              extent={{-31,-5},{28,-64}},
              lineColor={0,0,0},
              textString="dn")}));
    end MolarFlowSensor;

    model MoleFractionSensor "Measure of mole fraction"
      extends Modelica.Icons.RoundSensor;
      extends Interfaces.PartialSubstanceSensor;

      Modelica.Blocks.Interfaces.RealOutput moleFraction(final unit="1")
      "Mole fraction of the substance"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={-100,0},
          rotation=180)));

    equation
      port_a.q = 0;

      moleFraction = x;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
              textString="x"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MoleFractionSensor;

    model ElectroChemicalPotentialSensor
    "Measure of electro-chemical potential"
      extends Modelica.Icons.RoundSensor;

      Modelica.Blocks.Interfaces.RealOutput u(final unit="J/mol")
      "Electro-chemical potential of the substance"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={-100,0},
          rotation=180)));

    Interfaces.SubstancePort_b port_a
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));
    equation

      port_a.u = u;

      port_a.q = 0;
      port_a.h_outflow = 0;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="u"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ElectroChemicalPotentialSensor;

    model MolalitySensor "Measure of molality of the substance"
      extends Modelica.Icons.RoundSensor;
      extends Interfaces.PartialSubstanceSensor;

    parameter Modelica.Units.SI.AmountOfSubstance
      AmountOfSolutionPer1kgOfSolvent=1
      "Amount of all particles in the solution per one kilogram of solvent";

       Modelica.Blocks.Interfaces.RealOutput molality(final unit="mol/kg")
      "Molality of the substance (amount of substance per mass of solvent)"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={-100,0},
          rotation=180)));

  protected
    constant Modelica.Units.SI.Mass KG=1;
    equation
      port_a.q = 0;

      x=molality*KG / AmountOfSolutionPer1kgOfSolvent;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="b"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MolalitySensor;

    model MolarConcentrationSensor "Measure of molarity of the substance"
      extends Modelica.Icons.RoundSensor;
      extends Interfaces.PartialSubstanceSensor;

    parameter Modelica.Units.SI.AmountOfSubstance
      AmountOfSolutionInOneLiter=1
      "Amount of all particles in one liter of the solution";

       Modelica.Blocks.Interfaces.RealOutput molarConcentration(final unit="mol/m3", displayUnit="mol/l")
      "Molarity of the substance (amount of substance in one liter of whole solution)"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={-100,0},
          rotation=180)));

  protected
    constant Modelica.Units.SI.Volume L=0.001;
    equation
      port_a.q = 0;

      x=molarConcentration*L / AmountOfSolutionInOneLiter;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="c"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MolarConcentrationSensor;

    model MassFractionSensor "Measure of mass fraction of the substance"
      extends Modelica.Icons.RoundSensor;
      extends Interfaces.PartialSubstanceSensor;

    parameter Modelica.Units.SI.AmountOfSubstance
      AmountOfSolutionInOneKilogram=1
      "Amount of all particles in one kilogram of the solution";

       Modelica.Blocks.Interfaces.RealOutput massFraction(final unit="kg/kg")
      "Mass fraction of the substance (mass of the substance per mass of the whole solution)"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={-100,0},
          rotation=180)));

    equation
      port_a.q = 0;

      x=(massFraction*stateOfMatter.specificAmountOfParticles(substanceData)) / AmountOfSolutionInOneKilogram;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="mx"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MassFractionSensor;

    model PartialPressureSensor
    "Measure of partial pressure of the substance in gaseous solution"
      extends Modelica.Icons.RoundSensor;
      extends Interfaces.PartialSubstanceSensor;

       Modelica.Blocks.Interfaces.RealOutput partialPressure(final unit="Pa")
      "Partial pressure of the substance in gaseous solution"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={-100,0},
          rotation=180)));

    equation
      port_a.q = 0;

      partialPressure = x*solution.p;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="p"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PartialPressureSensor;

    model DissociationCoefficient
    "Meassure dissociation coefficient (mole fraction based) for pure substances"
      extends Modelica.Icons.RectangularSensor;

      outer Modelica.Fluid.System system "System wide properties";

      parameter Boolean useTemperatureInput = false
      "=true, if temperature is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.Units.SI.Temperature T=system.T_ambient
      "Temperature if not useTemperatureInput" annotation (HideResult=
          true, Dialog(enable=not useTemperatureInput));

      Modelica.Blocks.Interfaces.RealInput temperature(start=
            T, final unit="K")=_temperature if useTemperatureInput
      "Temperature"
        annotation (HideResult=true,Placement(transformation(extent={{-120,58},
                {-80,98}}), iconTransformation(extent={{-20,-20},{20,20}},
            rotation=270,
            origin={-60,40})));


      parameter Boolean useTotalAmountOfSubstancesInput = false
      "=true, if total amount of substances in solution is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.AmountOfSubstance n=1
      "Amount of all substances in solution per one liter of solution if not useTotalAmountOfSubstancesInput"
      annotation (HideResult=true, Dialog(enable=not
            useTotalAmountOfSubstancesInput));

      Modelica.Blocks.Interfaces.RealInput totalAmountOfSubstances(start=
            n, final unit="mol")=_n if useTotalAmountOfSubstancesInput
      "Temperature"
        annotation (HideResult=true,Placement(transformation(extent={{-120,58},
                {-80,98}}), iconTransformation(extent={{-20,-20},{20,20}},
            rotation=270,
            origin={40,40})));

    parameter Modelica.Units.SI.Mass m=1
      "Mass of solvent per one liter of solution";

      parameter Integer nS=0 "Number of substrates types"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, group="Ports"));

    parameter Modelica.Units.SI.StoichiometricNumber s[nS]=ones(nS)
      "Stoichiometric reaction coefficient for substrates"
      annotation (HideResult=true);

      parameter Integer nP=0 "Number of products types"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, group="Ports"));

    parameter Modelica.Units.SI.StoichiometricNumber p[nP]=ones(nP)
      "Stoichiometric reaction coefficients for products"
      annotation (HideResult=true);

    Interfaces.SubstancePort_b products[nP] "Products"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    Interfaces.SubstancePort_b substrates[nS] "Substrates"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

    Modelica.Units.SI.MolarEnergy DrG "Free Gibbs energy of reaction";

      Modelica.Blocks.Interfaces.RealOutput DissociationCoefficient_MoleFractionBased
      "Dissociation constant (if all substances has activity=1)"   annotation (Placement(transformation(
              extent={{-6,-86},{14,-66}}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-80})));

      Real DissociationCoefficient_MolalityBased
      "As ratio of molalities in moles per 1 kg of solvent";
      Real DissociationCoefficient_MolarityBased
      "As ratio of molar concentration in moles per liter of solution";

      Real pK
      "= -log10('mole-fraction based dissociation coefficient')";

  protected
    Modelica.Units.SI.Temperature _temperature;
    Modelica.Units.SI.AmountOfSubstance _n;
    equation
      if not useTemperatureInput then
        _temperature = T;
      end if;
      if not useTotalAmountOfSubstancesInput then
        _n = n;
      end if;

      substrates.q = zeros(nS);
      substrates.h_outflow = zeros(nS);

      products.q = zeros(nP);
      products.h_outflow = zeros(nP);


      DrG = ((p * products.u) - (s * substrates.u));

      DissociationCoefficient_MoleFractionBased = exp(-DrG/(Modelica.Constants.R*T));

      pK=-log10(DissociationCoefficient_MoleFractionBased);

      DissociationCoefficient_MolalityBased = ((n/m)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased;

      DissociationCoefficient_MolarityBased = ((n/1)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased;

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-160,-94},{-12,-68}},
              lineColor={0,0,0},
            textString="%s"),
            Text(
              extent={{12,-92},{160,-66}},
              lineColor={0,0,0},
            textString="%p")}),
        Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
<p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
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
    end DissociationCoefficient;

    model ActivityCoefficient
      "Calculate activity coefficient for product[1]"
      extends Modelica.Icons.RectangularSensor;

      outer Modelica.Fluid.System system "System wide properties";

      parameter Boolean useTemperatureInput = false
      "=true, if temperature is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Temperature T=system.T_ambient
      "Temperature if not useTemperatureInput" annotation (HideResult=
          true, Dialog(enable=not useTemperatureInput));

      Modelica.Blocks.Interfaces.RealInput temperature(start=
            T, final unit="K")=_temperature if useTemperatureInput
      "Temperature"
        annotation (HideResult=true,Placement(transformation(extent={{-120,58},
                {-80,98}}), iconTransformation(extent={{-20,-20},{20,20}},
            rotation=270,
            origin={-60,40})));

      parameter Boolean useTotalAmountOfSubstancesInput = false
      "=true, if total amount of substances in solution is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.AmountOfSubstance n=1
      "Amount of all substances in solution per one liter of solution if not useTotalAmountOfSubstancesInput"
      annotation (HideResult=true, Dialog(enable=not
            useTotalAmountOfSubstancesInput));

      Modelica.Blocks.Interfaces.RealInput totalAmountOfSubstances(start=
            n, final unit="mol")=_n if useTotalAmountOfSubstancesInput
      "Temperature"
        annotation (HideResult=true,Placement(transformation(extent={{-120,58},
                {-80,98}}), iconTransformation(extent={{-20,-20},{20,20}},
            rotation=270,
            origin={40,40})));

    parameter Modelica.Units.SI.Mass m=1
      "Mass of solvent per one liter of solution";

      parameter Integer nS=0 "Number of substrates types"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, group="Ports"));

    parameter Modelica.Units.SI.StoichiometricNumber s[nS]=ones(nS)
      "Stoichiometric reaction coefficient for substrates"
      annotation (HideResult=true);

      parameter Integer nP=0 "Number of products types"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, group="Ports"));

    parameter Modelica.Units.SI.StoichiometricNumber p[nP]=ones(nP)
      "Stoichiometric reaction coefficients for products"
      annotation (HideResult=true);

    Interfaces.SubstancePort_b products[nP] "Products"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    Interfaces.SubstancePort_b substrates[nS] "Substrates"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

    Modelica.Units.SI.MolarEnergy DrG "Free Gibbs energy of reaction";

      Modelica.Blocks.Interfaces.RealOutput activityCoeficient
      "Activity coeficient of one product"   annotation (Placement(transformation(
              extent={{-6,-86},{14,-66}}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-80})));

      parameter Boolean MolarityBased = true "if dissociation coefficient is molarity based";

      parameter Real DissociationCoefficient_MoleFractionBased = if MolarityBased then DissociationCoefficient_MolarityBased/((n/1)^(p*ones(nP)-s*ones(nS))) else DissociationCoefficient_MolalityBased/((n/m)^(p*ones(nP)-s*ones(nS)))
      "K as ratio of mole fractions";
      parameter Real DissociationCoefficient_MolalityBased = ((n/m)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased
      "K as ratio of molalities in moles per 1 kg of solvent"
      annotation (HideResult=true, Dialog(enable=not MolarityBased));
      parameter Real DissociationCoefficient_MolarityBased = ((n/1)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased
      "K as ratio of molar concentration in moles per liter of solution"
      annotation (HideResult=true, Dialog(enable=MolarityBased));

      Real pK
      "= -log10('mole-fraction based dissociation coefficient')";

  protected
    Modelica.Units.SI.Temperature _temperature;
    Modelica.Units.SI.AmountOfSubstance _n;
    equation
      if not useTemperatureInput then
        _temperature = T;
      end if;
      if not useTotalAmountOfSubstancesInput then
        _n = n;
      end if;

      substrates.q = zeros(nS);
      substrates.h_outflow = zeros(nS);

      products.q = zeros(nP);
      products.h_outflow = zeros(nP);

      DrG = ((p * products.u) - (s * substrates.u)) + (if (nP>0) then p[1] else 1)*Modelica.Constants.R*T*log(activityCoeficient);

      DissociationCoefficient_MoleFractionBased = exp(-DrG/(Modelica.Constants.R*T));

      pK=-log10(DissociationCoefficient_MoleFractionBased);

      //DissociationCoefficient_MolalityBased = ((n/m)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased;

      //DissociationCoefficient_MolarityBased = ((n/1)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased;

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-160,-94},{-12,-68}},
              lineColor={0,0,0},
            textString="%s"),
            Text(
              extent={{12,-92},{160,-66}},
              lineColor={0,0,0},
            textString="%p")}),
        Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
<p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
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
    end ActivityCoefficient;
  end Sensors;

  package Sources "Chemical sources"
    extends Modelica.Icons.SourcesPackage;

    model PureSubstance "Constant source of pure substance"
      extends Interfaces.PartialSubstance;

      parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
        "Temperature"
        annotation (HideResult=true, Dialog(enable=not useTemperatureInput));
      parameter Modelica.Units.SI.Pressure Pressure=system.p_ambient
        "Pressure"
        annotation (HideResult=true, Dialog(enable=not usePressureInput));
      parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
        "Electric potential" annotation (HideResult=true, Dialog(enable=not
              useElectricPotentialInput));
      parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
         0 "Ionic strength" annotation (HideResult=true, Dialog(enable=not
              useIonicStrengthInput));

      parameter Boolean useTemperatureInput = false
      "=true, if temperature is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Boolean usePressureInput = false
      "=true, if pressure is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Boolean useElectricPotentialInput = false
      "=true, if electric potential is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Boolean useIonicStrengthInput = false
      "=true, if ionic strength is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      Modelica.Blocks.Interfaces.RealInput T(start=
            Temperature, final unit="K")=temperature if useTemperatureInput
      "Temperature"
        annotation (HideResult=true,Placement(transformation(extent={{-120,58},
                {-80,98}}), iconTransformation(extent={{-120,58},{-80,98}})));

      Modelica.Blocks.Interfaces.RealInput p(start=
            Pressure, final unit="Pa")=pressure if usePressureInput
      "Pressure"
        annotation (HideResult=true,Placement(transformation(extent={{-120,16},
                {-80,56}}), iconTransformation(extent={{-120,16},{-80,56}})));

      Modelica.Blocks.Interfaces.RealInput v(start=
            ElectricPotential, final unit="Pa")=electricPotential if useElectricPotentialInput
      "Electric potential"
        annotation (HideResult=true,Placement(transformation(extent={{-120,-60},
                {-80,-20}}),iconTransformation(extent={{-120,-60},{-80,-20}})));

      Modelica.Blocks.Interfaces.RealInput I(start=
            MoleFractionBasedIonicStrength, final unit="mol/mol")=moleFractionBasedIonicStrength if useIonicStrengthInput
      "Pressure"
        annotation (HideResult=true,Placement(transformation(extent={{-120,-100},
                {-80,-60}}),iconTransformation(extent={{-120,-100},{-80,-60}})));
  protected
      Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
          Modelica.Constants.R*temperature))
        "Dissociation constant of hydrogen bond between base molecules";
      Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
          stateOfMatter.selfClusteringBondEnthalpy(
                                               substanceData) - temperature*
          stateOfMatter.selfClusteringBondEntropy(
                                              substanceData)
        "Gibbs energy of hydrogen bond between H2O molecules";


    equation


       if stateOfMatter.selfClustering(substanceData) then

        //Liquid cluster theory - equilibrium:
        //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules

        //sum(x[i]) = x/(1-K*x) = amountOfParticles/amountOfParticles = 1;
        x = 1/(1+SelfClustering_K) "mole fraction of free base molecule";
      else
        x = 1 "pure substance is composed only with free base molecules";
      end if;



       if (not useTemperatureInput) then
         temperature = Temperature;
       end if;
       if (not usePressureInput) then
         pressure = Pressure;
       end if;
       if (not useElectricPotentialInput) then
         electricPotential = ElectricPotential;
       end if;
       if (not useIonicStrengthInput) then
         moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;
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
              lineColor={128,0,255}),
            Text(
              extent={{-104,-76},{100,-100}},
              lineColor={0,0,0},
              textString="%T K")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PureSubstance;

    model ExternalIdealGasSubstance
    "Ideal gas substance with defined partial pressure"
      extends Interfaces.PartialSubstance(redeclare package stateOfMatter =
            Interfaces.IdealGas);

      parameter Boolean usePartialPressureInput = false
      "=true, if fixed partial pressure is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.Units.SI.Pressure PartialPressure=0
        "Fixed partial pressure if usePartialPressureInput=false" annotation (
         HideResult=true, Dialog(enable=not usePartialPressureInput));

      parameter Modelica.Units.SI.Pressure TotalPressure=system.p_ambient
        "Total pressure of the whole gaseous solution";

      parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
        "Temperature";
      parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
         0 "Ionic strength";
      parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput partialPressure(start=
            PartialPressure, final unit="Pa")=p if usePartialPressureInput
      "Partial pressure of gas = total pressure * gas fraction"
        annotation (HideResult=true,Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.Units.SI.Pressure p "Current partial pressure";

      parameter Modelica.Units.SI.Volume Volume=0.001
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
              lineColor={128,0,255}),
            Text(
              extent={{-100,-102},{104,-126}},
              lineColor={0,0,0},
              textString="%T K")}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ExternalIdealGasSubstance;

    model ExternalMolality "Constant source of substance molality"
      extends Interfaces.PartialSubstance;

      outer Modelica.Fluid.System system "System wide properties";

       parameter Real Molality(final unit="mol/kg") = 1e-8
      "Fixed molality of the substance if useMolalityInput=false"
        annotation (HideResult=true, Dialog(enable=not useMolalityInput));

      parameter Modelica.Units.SI.AmountOfSubstance
        AmountOfSolutionPer1KgSolvent=55.508
        "Amount of all particles in the solution per one kilogram of solvent";

        parameter Boolean useMolalityInput = false
      "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
        "Temperature";
      parameter Modelica.Units.SI.Pressure Pressure=system.p_ambient
        "Pressure";
      parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
         0 "Ionic strength";
      parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
        "Electric potential";



      Modelica.Blocks.Interfaces.RealInput molalityInput(start=Molality,final unit="mol/kg")=n/KG if
           useMolalityInput
        annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.Units.SI.AmountOfSubstance n "Current amount of the substance";

  protected
      constant Modelica.Units.SI.Mass KG=1;
    equation
       if not useMolalityInput then
         n=Molality*KG;
       end if;

      x = n/AmountOfSolutionPer1KgSolvent;

      //solution properties at the port
      temperature = Temperature;
      pressure = Pressure;
      electricPotential = ElectricPotential;
      moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;


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
              lineColor={128,0,255}),
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
    end ExternalMolality;

    model ExternalConcentration "Constant source of molar concentration"
       extends Interfaces.PartialSubstance;

       outer Modelica.Fluid.System system "System wide properties";

       parameter Real MolarConcentration(final unit="mol/m3", displayUnit="mol/l") = 1e-8
      "Fixed molarity of the substance if useMolarityInput=false"
        annotation (HideResult=true, Dialog(enable=not useMolarityInput));

      parameter Modelica.Units.SI.AmountOfSubstance AmountOfSolutionIn1L=55.508
        "Amount of all particles in the solution one liter of solvent";

        parameter Boolean useMolarityInput = false
      "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
        "Temperature";
      parameter Modelica.Units.SI.Pressure Pressure=system.p_ambient
        "Pressure";
      parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
         0 "Ionic strength";
      parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
        "Electric potential";

      Modelica.Blocks.Interfaces.RealInput molarConcentrationInput(start=MolarConcentration,final unit="mol/m3", displayUnit="mol/l")=n/L if
           useMolarityInput
        annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.Units.SI.AmountOfSubstance n "Current amount of the substance";

  protected
      constant Modelica.Units.SI.Volume L=0.001;
    equation
       if not useMolarityInput then
         n=MolarConcentration*L;
       end if;

      x = n/AmountOfSolutionIn1L;

      //solution properties at the port
      temperature = Temperature;
      pressure = Pressure;
      electricPotential = ElectricPotential;
      moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;


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
    end ExternalConcentration;

    model ExternalMoleFraction "Constant source of substance mole fraction"
         extends Interfaces.PartialSubstance;

       outer Modelica.Fluid.System system "System wide properties";

      parameter Modelica.Units.SI.MoleFraction MoleFraction=1e-8
        "Fixed mole fraction of the substance if useMoleFractionInput=false"
        annotation (HideResult=true, Dialog(enable=not useMoleFractionInput));

        parameter Boolean useMoleFractionInput = false
      "Is mole fraction of the substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
        "Temperature";
      parameter Modelica.Units.SI.Pressure Pressure=system.p_ambient
        "Pressure";
      parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
         0 "Ionic strength";
      parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
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
              lineColor={128,0,255}),
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
    end ExternalMoleFraction;

    model ExternalElectroChemicalPotential
    "Constant source of electro-chemical potential"

    parameter Modelica.Units.SI.ChemicalPotential U=1e-8
      "Fixed electro-chemical potential of the substance if usePotentialInput=false"
      annotation (HideResult=true, Dialog(enable=not usePotentialInput));

       parameter Boolean usePotentialInput = false
      "Is electro-chemical potential of the substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      Modelica.Blocks.Interfaces.RealInput uInput(final unit="J/mol")=port_a.u if
           usePotentialInput annotation (HideResult=true, Placement(transformation(
              extent={{-120,-20},{-80,20}})));

    Interfaces.SubstancePort_a port_a
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));
    parameter Modelica.Units.SI.MolarEnthalpy MolarHeat=0;
    equation
       if not usePotentialInput then
         port_a.u=U;
       end if;


      port_a.h_outflow = MolarHeat;


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
              lineColor={128,0,255}),
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
    end ExternalElectroChemicalPotential;

    model SubstanceInflow "Molar pump of substance to system"
      extends Interfaces.ConditionalSubstanceFlow;

    Interfaces.SubstancePort_b port_b "Outflow"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      parameter Modelica.Units.SI.MolarEnthalpy MolarHeat=0;
    equation
      port_b.q = -q;
      port_b.h_outflow = MolarHeat;

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
              extent={{-220,-80},{220,-60}},
              textString="%name",
              lineColor={128,0,255})}),      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceInflow;

    model SubstanceOutflow "Molar pump of substance out of system"
      extends Interfaces.ConditionalSubstanceFlow;

    Interfaces.SubstancePort_b port_a "Inflow"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

    parameter Modelica.Units.SI.MolarEnthalpy MolarHeat=0;
    equation
      port_a.q = q;

      port_a.h_outflow = MolarHeat;

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
              extent={{-220,-80},{220,-60}},
              textString="%name",
              lineColor={128,0,255})}),      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceOutflow;

    model Clearance "Physiological Clearance"
     extends Interfaces.ConditionalSolutionFlow(final SolutionFlow=Clearance/K);
     extends Interfaces.PartialSubstanceSensor;

    parameter Modelica.Units.SI.VolumeFlowRate Clearance=0
      "Physiological clearance of the substance if useSolutionFlowInput=false"
      annotation (HideResult=true, Dialog(enable=not useSolutionFlowInput));

      parameter Real K(unit="1")=1
      "Coefficient such that Clearance = K*solutionFlow";

    Modelica.Units.SI.MolarFlowRate molarClearance
      "Current molar clearance";

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
              lineColor={128,0,255}),
            Text(
              extent={{-100,-30},{100,-50}},
              lineColor={0,0,0},
              textString="K=%K")}),        Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Clearance;

    model Degradation "Degradation of the substance"
      extends Interfaces.PartialSubstanceSensor;

    parameter Modelica.Units.SI.Time HalfTime
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
              lineColor={128,0,255}),
            Text(
              extent={{-100,54},{100,28}},
              lineColor={0,0,0},
              textString="t1/2 = %HalfTime s"),
            Polygon(
              points={{54,24},{54,-24},{44,-22},{44,22},{54,24}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{30,20},{30,-20},{20,-18},{20,18},{30,20}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{8,16},{8,-16},{-2,-14},{-2,14},{8,16}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-12,12},{-12,-12},{-22,-10},{-22,10},{-12,12}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-34,8},{-34,-8},{-44,-6},{-44,6},{-34,8}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-56,4},{-56,-4},{-66,-2},{-66,2},{-56,4}},
              lineColor={0,0,127},
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
<td>2009-2020</td>
</tr>
</table>
</html>"));
    end Degradation;

    model Buffer
    "Source of substance bounded to constant amount of buffer to reach linear dependence between concentration and electrochemical potential"
      extends Icons.Buffer;
           extends Interfaces.PartialSubstanceInSolution(
                     a(start = a_start));

    parameter Modelica.Units.SI.MoleFraction a_start=1e-7
      "Initial value of mole fraction of the buffered substance";

    parameter Modelica.Units.SI.AmountOfSubstance BufferValue=0.001
      "Fixed buffer value (slope between amount of buffered substance and -log10(activity)) if useBufferValueInput=false"
      annotation (HideResult=true, Dialog(enable=not useBufferValueInput));

         parameter Boolean useBufferValueInput = false
      "Is buffer value of the substance an input?"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

          extends Interfaces.ConditionalKinetics(KC=1/(Modelica.Constants.R*298.15));

          Real bufferValue(final unit="1");

        Modelica.Blocks.Interfaces.RealInput bufferValueInput(
          final unit="mol/mol",
          start=BufferValue)=bufferValue if
             useBufferValueInput annotation (HideResult=true, Placement(transformation(
                extent={{-120,-20},{-80,20}})));

          Real xref;
    Modelica.Units.SI.AmountOfSubstance nFreeBuffer(start=-log10(a_start)
          *BufferValue) "amount of base molecules without H+";
    Modelica.Units.SI.MoleFraction xFreeBuffer;

  protected
    Modelica.Units.SI.MolarEnthalpy streamEnthalpy;

        constant Real InvLog_10=1/log(10);
    initial equation
        xFreeBuffer = -log10(a_start)*(bufferValue/solution.n);

    equation
        if not useBufferValueInput then
          bufferValue = BufferValue;
        end if;

        der(nFreeBuffer) = -port_a.q;
        // <- This is mathematically the same as two following lines. However, the differential solvers can handle the log10n much better. :-)
        //der(log10nFreeBuffer)=(InvLog_10)*(port_a.q/nFreeBuffer);
        //nFreeBuffer = 10^log10nFreeBuffer;

        xFreeBuffer = nFreeBuffer/solution.n;
       // port_a.q = (solution.n*KC)*(xFreeBuffer - xref);
        port_a.q = KC*(Modelica.Constants.R*solution.T*log(xFreeBuffer) - Modelica.Constants.R*solution.T*log(xref)); //alternative kinetics
        xref = -log10(a)*(bufferValue/solution.n);

      //solution flows
      streamEnthalpy = actualStream(port_a.h_outflow);

      solution.dH =streamEnthalpy*port_a.q - der(molarEnthalpy)*nFreeBuffer;
      solution.i = Modelica.Constants.F * z * port_a.q - Modelica.Constants.F*der(z)*nFreeBuffer;
      solution.dV = molarVolume * port_a.q - der(molarVolume)*nFreeBuffer;

      //extensive properties
      solution.nj=0;
      solution.mj=-nFreeBuffer*stateOfMatter.molarMassOfBaseMolecule(substanceData);
      solution.Vj=-nFreeBuffer*molarVolume;
      solution.Gj=-nFreeBuffer*port_a.u;
      solution.Qj=-Modelica.Constants.F*nFreeBuffer*z;
      solution.Ij=-(1/2) * ( nFreeBuffer * z^2);


        annotation ( Icon(coordinateSystem(
                preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
              graphics={
              Text(
                extent={{-82,62},{92,24}},
                textString="%name",
                lineColor={128,0,255})}),
          Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Buffer;

    model SubstanceInflowT
      "Molar pump of substance at defined temperature to system"
      extends Interfaces.ConditionalSubstanceFlow;

    Interfaces.SubstancePort_b port_b "Outflow"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));

     outer Modelica.Fluid.System system "System wide properties";

     replaceable package stateOfMatter =
        Chemical.Interfaces.Incompressible                                  constrainedby
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

     parameter stateOfMatter.SubstanceData substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true);

     parameter Modelica.Units.SI.Temperature T=system.T_ambient;
    equation
      port_b.q = -q;
      port_b.h_outflow = stateOfMatter.molarEnthalpy(substanceData,T);

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
              extent={{-220,-80},{220,-60}},
              textString="%name",
              lineColor={128,0,255})}),      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceInflowT;
  end Sources;

  package Interfaces "Chemical interfaces"
    extends Modelica.Icons.InterfacesPackage;

    connector SubstancePort
    "Electro-chemical potential and molar change of the substance in the solution"

    Modelica.Units.SI.ChemicalPotential u
      "Electro-chemical potential of the substance in the solution";

    flow Modelica.Units.SI.MolarFlowRate q
      "Molar change of the substance";

      //with molar flow of substance heat energy is changing also..
    stream Modelica.Units.SI.MolarEnthalpy h_outflow
      "Outgoing molar enthalphy";

      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
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
    end SubstancePort;

    connector SubstancePort_a
    "Electro-chemical potential and molar flow of the substance in the solution"
      extends SubstancePort;

    annotation (
        defaultComponentName="port_a",
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                100}}),     graphics={Rectangle(
              extent={{-20,10},{20,-10}},
              lineColor={158,66,200}),Rectangle(
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
    end SubstancePort_a;

    connector SubstancePort_b
    "Electro-chemical potential and molar flow of the substance in the solution"
      extends SubstancePort;

    annotation (
        defaultComponentName="port_b",
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                100}}),     graphics={Rectangle(
              extent={{-20,10},{20,-10}},
              lineColor={158,66,200}),Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={158,66,200},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}),
            graphics={Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={158,66,200},
              lineThickness=1,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
       Text(extent = {{-160,110},{40,50}}, lineColor={172,72,218},   textString = "%name")}),
        Documentation(info="<html>
<p>Chemical port with external definition of the substance outside the component.</p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstancePort_b;

    partial model OnePort "Base model for chemical process"

    SubstancePort_b port_a annotation (Placement(transformation(extent={{-110,-10},
              {-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));
    SubstancePort_a port_b annotation (Placement(transformation(extent={{90,-10},
              {110,10}}), iconTransformation(extent={{90,-10},{110,10}})));

    parameter Boolean EnthalpyNotUsed=false annotation (
        Evaluate=true,
        HideResult=true,
        choices(checkBox=true),
        Dialog(tab="Advanced", group="Performance"));

    equation
      port_a.q + port_b.q = 0;
      port_a.h_outflow =
       if EnthalpyNotUsed then 0
       else inStream(port_b.h_outflow);
      port_b.h_outflow =
       if EnthalpyNotUsed then 0
       else inStream(port_a.h_outflow);
    end OnePort;

    connector SubstancePorts_a
      extends SubstancePort;
      annotation (
         defaultComponentName="ports_a",
         Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-50,-200},{50,200}},
            initialScale=0.2),graphics={
            Text(extent={{-73,130},{77,100}}, textString="%name"),
            Rectangle(
              extent={{25,-100},{-25,100}},
              lineColor={158,66,200}),
                      Rectangle(
              extent={{-20,20},{20,-20}},
              lineColor={158,66,200},
              lineThickness=1),
                      Rectangle(
              extent={{-20,90},{20,50}},
              lineColor={158,66,200},
              lineThickness=1),
                      Rectangle(
              extent={{-20,-52},{20,-90}},
              lineColor={158,66,200},
              lineThickness=1)}),
               Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-50,-200},{50,200}},
            initialScale=0.2),graphics={
            Rectangle(
              extent={{50,-200},{-50,200}},
              lineColor={158,66,200},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
                                      Rectangle(
              extent={{-40,38},{40,-42}},
              lineColor={158,66,200},
              fillColor={158,66,200},
              fillPattern=FillPattern.Solid),
                                      Rectangle(
              extent={{-40,170},{40,90}},
              lineColor={158,66,200},
              fillColor={158,66,200},
              fillPattern=FillPattern.Solid),
                                      Rectangle(
              extent={{-40,-92},{40,-172}},
              lineColor={158,66,200},
              fillColor={158,66,200},
              fillPattern=FillPattern.Solid)}));

    end SubstancePorts_a;

    connector SubstancePorts_b
      extends SubstancePort;
      annotation (
         defaultComponentName="ports_b",
         Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-50,-200},{50,200}},
            initialScale=0.2),graphics={
            Text(extent={{-73,130},{77,100}}, textString="%name"),
            Rectangle(
              extent={{25,-100},{-25,100}},
              lineColor={158,66,200}),
                      Rectangle(
              extent={{-20,20},{20,-20}},
              lineColor={158,66,200},
              lineThickness=1),
                      Rectangle(
              extent={{-20,90},{20,50}},
              lineColor={158,66,200},
              lineThickness=1),
                      Rectangle(
              extent={{-20,-52},{20,-90}},
              lineColor={158,66,200},
              lineThickness=1)}),
               Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-50,-200},{50,200}},
            initialScale=0.2),graphics={
            Rectangle(
              extent={{50,-200},{-50,200}},
              lineColor={158,66,200},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
                                      Rectangle(
              extent={{-40,38},{40,-42}},
              lineColor={158,66,200}),Rectangle(
              extent={{-40,170},{40,90}},
              lineColor={158,66,200}),Rectangle(
              extent={{-40,-92},{40,-172}},
              lineColor={158,66,200})}));

    end SubstancePorts_b;

    partial model PartialSubstance

     outer Modelica.Fluid.System system "System wide properties";

      SubstancePort_a port_a "The substance"
     annotation (Placement(transformation(extent={{90,-10},{110,10}})));

     replaceable package stateOfMatter = Incompressible constrainedby StateOfMatter
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


     parameter stateOfMatter.SubstanceData substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true);

    Modelica.Units.SI.MoleFraction x "Mole fraction of the substance";

    Modelica.Units.SI.ActivityOfSolute a
      "Activity of the substance (mole-fraction based)";

  protected
    Modelica.Units.SI.ActivityCoefficient gamma
      "Activity coefficient of the substance";

    Modelica.Units.SI.ChargeNumberOfIon z "Charge number of ion";

    Modelica.Units.SI.Temperature temperature
      "Temperature of the solution";

    Modelica.Units.SI.Pressure pressure "Pressure of the solution";

    Modelica.Units.SI.ElectricPotential electricPotential
      "Electric potential of the solution";

    Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength
      "Ionic strength of the solution";

    //Modelica.Units.SI.MolarMass molarMass "Molar mass of the substance";

    Modelica.Units.SI.MolarEnthalpy molarEnthalpy
      "Molar enthalpy of the substance";

    Modelica.Units.SI.MolarEntropy molarEntropyPure
      "Molar entropy of the pure substance";

    Modelica.Units.SI.ChemicalPotential u0
      "Chemical potential of the pure substance";

    Modelica.Units.SI.ChemicalPotential uPure
      "Electro-Chemical potential of the pure substance";

    Modelica.Units.SI.MolarVolume molarVolume
      "Molar volume of the substance";

    Modelica.Units.SI.MolarVolume molarVolumePure
      "Molar volume of the pure substance";

    Modelica.Units.SI.MolarVolume molarVolumeExcess
      "Molar volume excess of the substance in solution (typically it is negative as can be negative)";

      //  Modelica.SIunits.MolarHeatCapacity molarHeatCapacityCp
      //    "Molar heat capacity of the substance at constant pressure";



    equation
     //aliases
     gamma = stateOfMatter.activityCoefficient(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     z = stateOfMatter.chargeNumberOfIon(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
    // molarMass = stateOfMatter.molarMass(substanceData);

     molarEnthalpy = stateOfMatter.molarEnthalpy(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     molarEntropyPure = stateOfMatter.molarEntropyPure(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     u0 = stateOfMatter.chemicalPotentialPure(
       substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength);
     uPure = stateOfMatter.electroChemicalPotentialPure(
       substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength);
     molarVolume = stateOfMatter.molarVolume(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     molarVolumePure = stateOfMatter.molarVolumePure(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     molarVolumeExcess = stateOfMatter.molarVolumeExcess(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     //  molarHeatCapacityCp = stateOfMatter.molarHeatCapacityCp(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);

     //activity of the substance
     a = gamma*x;

     //electro-chemical potential of the substance in the solution
     port_a.u = stateOfMatter.chemicalPotentialPure(
       substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength)
       + Modelica.Constants.R*temperature*log(a)
       + z*Modelica.Constants.F*electricPotential;

     port_a.h_outflow = molarEnthalpy;


     annotation (
       Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PartialSubstance;

    partial model PartialSubstanceInSolution
    "Substance properties for components, where the substance is connected with the solution"

      SolutionPort            solution
      "To connect substance with solution, where is pressented"                                  annotation (Placement(transformation(
              extent={{-70,-110},{-50,-90}}),iconTransformation(extent={{-70,-110},{
                -50,-90}})));

      extends PartialSubstance;

  protected
    Modelica.Units.SI.AmountOfSubstance amountOfSolution
      "Amount of all solution particles";

    equation

      //aliases
      temperature = solution.T;
      pressure = solution.p;
      electricPotential = solution.v;
      amountOfSolution = solution.n;
      moleFractionBasedIonicStrength = solution.I;


    end PartialSubstanceInSolution;

    partial model PartialSubstanceInSolutionWithAdditionalPorts
      "Substance properties for components, where the substance is connected with the solution"

      extends PartialSubstanceInSolution;

    Modelica.Units.SI.MolarFlowRate q
      "Molar flow rate of the substance into the component";

      SubstanceMassPort_a   port_m "Substance mass fraction port"
            annotation (Placement(transformation(extent={{92,-110},{112,-90}})));

      SubstanceMolarityPort_a
                            port_c
        annotation (Placement(transformation(extent={{90,90},{110,110}})));

    equation
      //molar mass flow
      q=(port_a.q + port_c.q + port_m.m_flow/stateOfMatter.molarMassOfBaseMolecule(substanceData));

      //substance mass fraction
      port_m.x_mass = solution.mj/solution.m;
      port_c.c = solution.nj/solution.V;

    end PartialSubstanceInSolutionWithAdditionalPorts;

    partial model PartialSubstanceSensor
    "Base class for sensor based on substance and solution properties"
      extends PartialSubstanceInSolution;

    equation
      //solution is not changed by the sensor components
      solution.dH = 0;
      solution.i = 0;
      solution.dV = 0;
      solution.Gj = 0;
      solution.nj = 0;
      solution.mj = 0;
      solution.Qj = 0;
      solution.Ij = 0;
      solution.Vj = 0;

    end PartialSubstanceSensor;

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
        output Real specificAmountOfSubstance(unit="mol/kg") "Amount of substance particles per its mass";
    protected
        Modelica.Units.SI.MolarEnergy SelfClustering_dG;
        Real SelfClustering_K;
      algorithm
        if not selfClustering(substanceData) then
          specificAmountOfSubstance := 1/substanceData.MolarWeight;
        else
          SelfClustering_dG :=selfClusteringBondEnthalpy(substanceData) - T*
            selfClusteringBondEntropy(substanceData);

          SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));

          specificAmountOfSubstance := 1/((SelfClustering_K + 1)*substanceData.MolarWeight);
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
        output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
          "Amount of substance free base molecule per substance mass";
    protected
        Modelica.Units.SI.MolarEnergy SelfClustering_dG;
        Real SelfClustering_K;
      algorithm
        if not selfClustering(substanceData) then
          specificAmountOfFreeBaseMolecule := 1/substanceData.MolarWeight;
        else
          SelfClustering_dG :=selfClusteringBondEnthalpy(substanceData) - T*
            selfClusteringBondEntropy(substanceData);

          SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
          specificAmountOfFreeBaseMolecule :=(1- SelfClustering_K/(SelfClustering_K + 1)) / substanceData.MolarWeight;
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
</html>",     info="<html>
<p>Solution port integrates all substances of the solution:</p>
<p>Such as if there are connected together with electric port, thermal port and with port composed with the amont of substance and molar change of substance.</p>
</html>"), Icon(graphics={            Rectangle(
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
          Chemical.Interfaces.StateOfMatter
        constrainedby StateOfMatter
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
</html>", info="<html>
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
</html>", info="<html>
<h4>amountOfSolution = &sum; amountOfSubstances</h4>
<h4>mass = &sum; massOfSubstances</h4>
<h4>volume = &sum; volumeOfSubstances</h4>
<h4>freeGibbsEnergy = &sum; freeGibbsEnergiesOfSubstances</h4>
<p>To calculate the sum of extensive substance's properties is misused the Modelica \"flow\" prefix even there are not real physical flows. </p>
</html>"));
    end PartialSolutionWithHeatPort;

    partial model ConditionalSolutionFlow
    "Input of solution molar flow vs. parametric solution molar flow"

      parameter Boolean useSolutionFlowInput = false
      "=true, if solution flow is provided via input"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),
              Dialog(group="Conditional inputs", __Dymola_compact=true));

    parameter Modelica.Units.SI.VolumeFlowRate SolutionFlow=0
      "Volume flow rate of the solution if useSolutionFlowInput=false"
      annotation (HideResult=true, Dialog(enable=not useSolutionFlowInput));

    parameter Modelica.Units.SI.AmountOfSubstance AmountOfSolutionIn1L=
        55.508 "The amount of all particles in one liter of the solution";

      Modelica.Blocks.Interfaces.RealInput solutionFlow(start=SolutionFlow, final unit="m3/s")=
         q*OneLiter/AmountOfSolutionIn1L if useSolutionFlowInput
         annotation ( HideResult=true, Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,40}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,40})));

    Modelica.Units.SI.MolarFlowRate q "Current molar solution flow";

  protected
    constant Modelica.Units.SI.Volume OneLiter=0.001 "One liter";

    equation
      if not useSolutionFlowInput then
        q*OneLiter/AmountOfSolutionIn1L = SolutionFlow;
      end if;

    end ConditionalSolutionFlow;

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

      Modelica.Blocks.Interfaces.RealInput substanceFlow(start=SubstanceFlow, final unit="mol/s")=q if
           useSubstanceFlowInput
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

      parameter Boolean useKineticsInput = false
        "= true, if kinetics coefficient is provided via input"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),
          Dialog(group="Chemical kinetics", __Dymola_compact=true));

      parameter Real KC(final unit="mol2.s-1.J-1")=1
        "Chemical kinetics coefficient if useKineticsInput=false"
        annotation (HideResult=true, Dialog(group="Chemical kinetics", enable=not useKineticsInput));

      Modelica.Blocks.Interfaces.RealInput kineticsCoefficientInput(start=KC, final unit="mol2.s-1.J-1")=
         kC if useKineticsInput
         annotation ( HideResult=true, Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={-60,40}),
                            iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={-60,40})));

      Real kC(final unit="mol2.s-1.J-1") "Current kinetics coefficient";

    equation
      if not useKineticsInput then
        kC = KC;
      end if;

    end ConditionalKinetics;

   connector SubstanceMassPort

    Modelica.Units.SI.MassFraction x_mass
      "Mass fraction of the substance in the solution";

    flow Modelica.Units.SI.MassFlowRate m_flow
      "Mass flow rate of the substance";
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
   end SubstanceMassPort;

     connector SubstanceMassPort_a
        "Mass fraction and mass flow of the substance in the solution"
        extends SubstanceMassPort;

      annotation (
          defaultComponentName="port_a",
          Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                  100}}),     graphics={Rectangle(
                extent={{-20,10},{20,-10}},
                lineColor={105,44,133}),Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={105,44,133},
                fillColor={105,44,133},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}),
              graphics={Rectangle(
                extent={{-40,40},{40,-40}},
                lineColor={105,44,133},
                fillColor={105,44,133},
                fillPattern=FillPattern.Solid,
                lineThickness=1),
         Text(extent = {{-160,110},{40,50}}, lineColor={105,44,133},   textString = "%name")}),
          Documentation(info="<html>
<p>Chemical port with internal definition of the substance inside the component. </p>
</html>", revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
     end SubstanceMassPort_a;

    connector SubstanceMassPort_b
      "Mass fraction and mass flow of the substance in the solution"
      extends SubstanceMassPort;

    annotation (
        defaultComponentName="port_b",
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                100}}),     graphics={Rectangle(
              extent={{-20,10},{20,-10}},
              lineColor={105,44,133}),Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={105,44,133},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}),
            graphics={Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={105,44,133},
              lineThickness=1,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
       Text(extent = {{-160,110},{40,50}}, lineColor={105,44,133},   textString = "%name")}),
        Documentation(info="<html>
<p>Chemical port with external definition of the substance outside the component.</p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceMassPort_b;

    connector SubstanceMassPorts_a
      extends SubstanceMassPort;
      annotation (
         defaultComponentName="ports_a",
         Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-50,-200},{50,200}},
            initialScale=0.2),graphics={
            Text(extent={{-73,130},{77,100}},
              textString="%name",
              lineColor={105,44,133}),
            Rectangle(
              extent={{25,-100},{-25,100}},
              lineColor={105,44,133}),
                      Rectangle(
              extent={{-20,20},{20,-20}},
              lineColor={105,44,133},
              lineThickness=1),
                      Rectangle(
              extent={{-20,90},{20,50}},
              lineColor={105,44,133},
              lineThickness=1),
                      Rectangle(
              extent={{-20,-52},{20,-90}},
              lineColor={105,44,133},
              lineThickness=1)}),
               Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-50,-200},{50,200}},
            initialScale=0.2),graphics={
            Rectangle(
              extent={{50,-200},{-50,200}},
              lineColor={105,44,133},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
                                      Rectangle(
              extent={{-40,38},{40,-42}},
              lineColor={105,44,133},
              fillColor={105,44,133},
              fillPattern=FillPattern.Solid),
                                      Rectangle(
              extent={{-40,170},{40,90}},
              lineColor={105,44,133},
              fillColor={105,44,133},
              fillPattern=FillPattern.Solid),
                                      Rectangle(
              extent={{-40,-92},{40,-172}},
              lineColor={105,44,133},
              fillColor={105,44,133},
              fillPattern=FillPattern.Solid)}));

    end SubstanceMassPorts_a;

    connector SubstanceMolarityPort

    Modelica.Units.SI.Concentration c
      "Molarity of the substance in the solution";

    flow Modelica.Units.SI.MolarFlowRate q
      "Molar flow rate of the substance";
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end SubstanceMolarityPort;

    connector SubstanceMolarityPort_a
      "Electro-chemical potential and molar flow of the substance in the solution"
      extends SubstanceMolarityPort;

    annotation (
        defaultComponentName="port_a",
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                100}}),     graphics={Rectangle(
              extent={{-20,10},{20,-10}},
              lineColor={174,73,220}),Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={174,73,220},
              fillColor={174,73,220},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}),
            graphics={Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={174,73,220},
              fillColor={174,73,220},
              fillPattern=FillPattern.Solid,
              lineThickness=1),
       Text(extent = {{-160,110},{40,50}}, lineColor={174,73,220},   textString = "%name")}),
        Documentation(info="<html>
<p>Chemical port with internal definition of the substance inside the component. </p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceMolarityPort_a;

    connector SubstanceMolarityPort_b
      "Electro-chemical potential and molar flow of the substance in the solution"
      extends SubstanceMolarityPort;

    annotation (
        defaultComponentName="port_b",
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
                100}}),     graphics={Rectangle(
              extent={{-20,10},{20,-10}},
              lineColor={174,73,220}),Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={174,73,220},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}),
            graphics={Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={174,73,220},
              lineThickness=1,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
       Text(extent = {{-160,110},{40,50}}, lineColor={174,73,220},   textString = "%name")}),
        Documentation(info="<html>
<p>Chemical port with external definition of the substance outside the component.</p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceMolarityPort_b;
  end Interfaces;

  annotation (
preferredView="info",
version="1.4.0",
versionDate="2021-01-27",
dateModified = "2021-01-27 11:10:41Z",
conversion(
  from(version="1.3.1", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.3_to_1.4.mos",
        to="1.4.0"),
  from(version="1.3.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.3_to_1.4.mos",
        to="1.4.0"),
  from(version="1.2.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.3_to_1.4.mos",
        to="1.4.0"),
  from(version="1.1.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.1_to_1.4.mos",
        to="1.4.0"),
  from(version="1.0.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.0_to_1.4.mos",
        to="1.4.0")),
      uses( Modelica(version="4.0.0")),
  Documentation(revisions="<html>
<p>Copyright (c) 2021, Marek Matej&aacute;k, Ph.D. </p>
<p>All rights reserved. </p>
<p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: </p>
<ol>
<li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. </li>
<li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. </li>
<li>Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. </li>
</ol>
<p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS &quot;AS IS&quot; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
</html>", info="<html>
<p>During each electro-chemical process an <a href=\"modelica://Chemical.Components.Substance\">electro-chemical potential</a> of the substances is equilibrating and all thermodynamical properties of the homogenous chemical solutions are evaluated. </p>
<p>Processes: chemical reactions, gas dissolution, diffusion, membrane transports, osmotic fluxes, electrochemical cells, electrodes, ..</p>
<p>Please see the <a href=\"modelica://Chemical.UsersGuide.Overview\">overview</a>.</p>
</html>"));
end Chemical;
