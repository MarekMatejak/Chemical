within Chemical;
package Obsolete
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
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.Solution\">Chemical solution</a></p><p>The solution is the base component of each model, because it defines the conditions of the electro-chemical processes. It integrates the total amount of substance (called amount of solution), heat, charge, entropy, volume and others from each substances to present the base properties such as temperature, pressure, electric potential and others. The usage is very simple - just connect each chemical substance with its chemical solution using their <a href=\"modelica://Chemical.Obsolete.Interfaces.SolutionPort\">SolutionPort</a>.</p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Substance1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.Substance\">Chemical substance</a></p><p>The chemical substance integrates the amount of the chemical substance and from the properties of the connected solution it presents the electro-chemical potential of the substance using the <a href=\"modelica://Chemical.Obsolete.Interfaces.SubstancePort\">SubstancePort</a>.</p><p>There are two basic <a href=\"modelica://Chemical.Obsolete.Interfaces.StateOfMatter\">states of matter</a>: <a href=\"modelica://Chemical.Obsolete.Interfaces.IdealGas\">ideal gas</a> and <a href=\"modelica://Chemical.Obsolete.Interfaces.Incompressible\">incompressible</a> substance. However, the user can easily (re)define their own state of matter by inserting the correct expressions for the pure substance <a href=\"modelica://Chemical.Obsolete.Interfaces.StateOfMatter.activityCoefficient\">activity coefficient</a>, <a href=\"modelica://Chemical.Obsolete.Interfaces.StateOfMatter.molarVolumePure\">molar volume</a>, <a href=\"modelica://Chemical.Obsolete.Interfaces.StateOfMatter.molarEntropyPure\">molar entropy</a> and <a href=\"modelica://Chemical.Obsolete.Interfaces.StateOfMatter.molarEnthalpyElectroneutral\">molar enthalpy</a>, based on the current solution state (temperature, pressure, electric potential and ionic strength) and the <a href=\"modelica://Chemical.Obsolete.Interfaces.StateOfMatter.SubstanceData\">substance data</a>. The object-oriented design allows users to define the substance data record as part of the state of matter package. Users can select substance parameters according to the state of matter, redefining the getter functions of substance properties.</p><p>The examples work with ideal gases in case of all gaseous substance and incompressible state of matter in case of liquid or solid. The definition data are the molar mass of the substance, the number of charges of the substance, the molar heat capacity of the substance at a constant pressure, free formation enthalpy, free formation Gibbs energy and density (if incompressible) &mdash; all at a temperature of 25&deg;C and pressure 1 bar. Since these parameters are usually recorded in chemical tables at this standard conditions. In this manner, more than 35 real chemical <a href=\"modelica://Chemical.Obsolete.Examples.Substances\">substances</a> in the example package of this chemical library have already been defined. The usage of these predefined substances&rsquo; data is very simple. In the parameter dialog of the chemical substance, the correct record with this data can be selected, as shown in Figure 1.</p><p>This setting is typically the most important setting of each chemical model. All equilibrium coefficients, standard voltages, dissolution coefficients, saturated vapor pressures and so on, are automatically solved using these substance data. As a result, for example, the chemical reaction component only needs to define the stoichiometry coefficients, and the connected substances reach equilibrium at the correct equilibrium coefficient.</p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Reaction1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.Reaction\">Chemical reaction</a></p><p>The chemical reaction component is very general. The dissociation constant of the equilibrium is calculated from substance properties at usual in thermodynamics, for example as definition of <a href=\"http://goldbook.iupac.org/S05915.html\">UIPAC</a>. For example if we want to define <a href=\"modelica://Chemical.Obsolete.Examples.SimpleReaction\">simple reaction A&lt;-&gt;B</a> with dissociation constant [B]/[A]=2 then it must be the difference between Gibbs energies of formation equal to B.DfG - A.DfG = - R * T * ln(2). Without lost of generality it is possible to select some substances as reference and give them the zero Gibbs energy of formation. The next substances created by some chemical process can be expressed from them such as example of <a href=\"modelica://Chemical.Obsolete.Examples.Hemoglobin.Allosteric_Hemoglobin_MWC\">alosteric hemoglobin</a> calculation. The kinetics of the chemical reaction is different as usual. However the most of processes can be recalculated with sufficient precision, for example the <a href=\"modelica://Chemical.Obsolete.Examples.EnzymeKinetics\">Michaelic-Menton</a> can be recalculated with precision of 1.5% of maximal rate. </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Diffusion1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.Diffusion\">Diffusion</a></p><p>Diffusion is a dynamic chemical process, wich is also equilibrating of electro-chemical potential of the substance. Analogically as in chemical reaction the speed of diffucion can be calculated as coefficient C multiplied by electro-chemical gratient. C can be a parammeter or input expressed from distance, substance and solution properties. </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/GasSolubility1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.GasSolubility\">Henry&apos;s law, Raoult&apos;s law or Sieverts&apos; law</a></p><p>Surprisingly, all these laws has the same basis = equilibrium of electro-chemical potential. The most of problems in data is caused by wrong selection of standard state as 1 mol/kg or 1 mol/L. Please avoid these assumptions of these totally confused states and use only mole fractions instead of each molality or molarity - the world will be much better (I promise). </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Membrane1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.Membrane\">Semipermeable membrane</a></p><p>The same as before - just equilibrating the electro-chemical potentials. A result is the Donnan&apos;s equilibrium, Nernst potentials of the ions and the membrane electric potential. Transporting water through membrane is reaching the osmotic equilibrium (The real one, not the simplified one defined by osmotic pressure lineary dependent on impermeable substance concentration). </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Speciation1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.Speciation\">Chemical speciation</a></p><p>The chemical speciation is for macromolecule composed with independent subunits is specific conformations. For example the hemoglobin is tetramer, which can be in two conformation: relaxed and tensed. In each of this conformation it has different afinities (different dissociation constant) for binding oxygen in each of four independent subunits. This alosteric effect can be modeled using speciation such as in <a href=\"modelica://Chemical.Obsolete.Examples.Hemoglobin.Allosteric_Hemoglobin2_MWC\">Allosteric_Hemoglobin2_MWC</a>. However the result should be the same as using the detailed reaction model <a href=\"modelica://Chemical.Obsolete.Examples.Hemoglobin.Allosteric_Hemoglobin_MWC\">Allosteric_Hemoglobin_MWC</a>.</p></td>
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
<td valign=\"middle\"><p><br><a href=\"modelica://Chemical.Obsolete.Interfaces.SubstancePort\">Chemical.Interfaces.SubstancePort</a> </p></td>
<td valign=\"middle\"><p><img src=\"modelica://Chemical/Resources/Images/UsersGuide/ChemicalPorts.png\"/></p></td>
</tr>
<tr>
<td valign=\"middle\"><h4>solution</h4></td>
<td valign=\"middle\"><p>p .. pressure of the solution</p><p>T .. temperature of the solution</p><p>v .. electric potential of the solution</p><p><br>n .. amount of all substances in the solution</p><p>m .. mass of the solution</p><p>V .. volume of the solution</p><p>G .. free Gibbs energy of the solution</p><p>Q .. electric charge of the solution</p><p>I .. ionic strength of the solution</p></td>
<td valign=\"middle\"><p>dV .. change of the volume of the solution</p><p>dH .. enthalpy change of the solution</p><p>i .. electric charge change of the solution</p><p><br><i>nj ..  amount of the substance</i></p><p><i>mj .. mass of the substance</i></p><p><i>Vj .. volume of the substance</i></p><p><i>Gj .. free Gibbs energy of the substance</i></p><p><i>Qj .. electric charge of the substance</i></p><p><i>Ij .. ionic strength of the substance</i></p></td>
<td valign=\"middle\"></td>
<td valign=\"middle\"><p><br><a href=\"modelica://Chemical.Obsolete.Interfaces.SolutionPort\">Chemical.Interfaces.SolutionPort</a></p></td>
<td valign=\"middle\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/SolutionPort.png\"/></p></td>
</tr>
<tr>
<td valign=\"middle\"><h4>substanceMass</h4></td>
<td valign=\"middle\"><p>x_mass .. mass fraction of the chemical substance in solution</p></td>
<td valign=\"middle\"><p>m_flow .. mass flow of the chemical substance</p></td>
<td valign=\"middle\"></td>
<td valign=\"middle\"><p><br><a href=\"modelica://Chemical.Obsolete.Interfaces.SubstanceMassPort\">Chemical.Interfaces.SubstanceMassPort</a> </p></td>
<td valign=\"middle\"><p><img src=\"modelica://Chemical/Resources/Images/UsersGuide/ChemicalMassPorts.png\"/></p></td>
</tr>
<tr>
<td valign=\"middle\"><h4>substanceMolarity</h4></td>
<td valign=\"middle\"><p>c .. molar concentration per liter of the chemical substance in solution</p></td>
<td valign=\"middle\"><p>q .. molar flow of the chemical substance</p></td>
<td valign=\"middle\"></td>
<td valign=\"middle\"><p><br><a href=\"modelica://Chemical.Obsolete.Interfaces.SubstanceMolarityPort\">Chemical.Interfaces.SubstanceMolarityPort</a> </p></td>
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

  class Version_1_4_1 "Version 1.4.1 (Nov. 23, 2023)"
    extends Modelica.Icons.ReleaseNotes;

  annotation (Documentation(info="<html>
<ul>
<li>FIX of Chemical.Interfaces.StateOfMatter.specificAmountOfFreeBaseMolecule,<br/>
  e.g. Chemical.Interfaces.Incompressible.specificAmountOfFreeBaseMolecule(Chemical.Substances.Water_liquid());
 = 0.018285267550297114</li>
</ul>
</html>"));
  end Version_1_4_1;
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
<p>5.května 50,</p>
<p>250 82 Úvaly u Prahy, </p>
<p>Czech Republic,</p>
<p>email: marek@matfyz.cz</p>
<h4><span style=\"color:#008000\">Organization: </span></h4>
<p>Institute of Pathological Physiology, First Faculty of Medicine, Charles University in Prague,</p>
<p>U Nemocnice 5, 128 53 Prague 2, Czech Republic</p>
<br><h4>Copyright notices of the files:</h4>
<p>Copyright (c) 2008-2023, Marek Matej&aacute;k, Charles University in Prague</p>
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

    model ElectronTransfer
    "Electron transfer from the solution to electric circuit"
      extends Icons.ElectronTransfer;
      extends Interfaces.PartialSubstanceInSolution(redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible, final substanceData=
            Chemical.Obsolete.Interfaces.Incompressible.SubstanceData(
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

      replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby Interfaces.StateOfMatter
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

      replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby Interfaces.StateOfMatter
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

  package Substances "Definitions of substances"
    package IdealGasesMSL

    record Ag "Ag(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ag);
     annotation (preferredView = "info");
    end Ag;

    record Agplus "Agplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Agplus,
      z=1);
     annotation (preferredView = "info");
    end Agplus;

    record Agminus "Agminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Agminus,
      z=-1);
     annotation (preferredView = "info");
    end Agminus;

    record Air "Air(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Air);
     annotation (preferredView = "info");
    end Air;

    record AL "AL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL);
     annotation (preferredView = "info");
    end AL;

    record ALplus "ALplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALplus,
      z=1);
     annotation (preferredView = "info");
    end ALplus;

    record ALminus "ALminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALminus,
      z=-1);
     annotation (preferredView = "info");
    end ALminus;

    record ALBr "ALBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALBr);
     annotation (preferredView = "info");
    end ALBr;

    record ALBr2 "ALBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALBr2);
     annotation (preferredView = "info");
    end ALBr2;

    record ALBr3 "ALBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALBr3);
     annotation (preferredView = "info");
    end ALBr3;

    record ALC "ALC(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALC);
     annotation (preferredView = "info");
    end ALC;

    record ALC2 "ALC2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALC2);
     annotation (preferredView = "info");
    end ALC2;

    record ALCL "ALCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALCL);
     annotation (preferredView = "info");
    end ALCL;

    record ALCLplus "ALCLplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALCLplus,
      z=1);
     annotation (preferredView = "info");
    end ALCLplus;

    record ALCL2 "ALCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALCL2);
     annotation (preferredView = "info");
    end ALCL2;

    record ALCL3 "ALCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALCL3);
     annotation (preferredView = "info");
    end ALCL3;

    record ALF "ALF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF);
     annotation (preferredView = "info");
    end ALF;

    record ALFplus "ALFplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALFplus,
      z=1);
     annotation (preferredView = "info");
    end ALFplus;

    record ALFCL "ALFCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALFCL);
     annotation (preferredView = "info");
    end ALFCL;

    record ALFCL2 "ALFCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALFCL2);
     annotation (preferredView = "info");
    end ALFCL2;

    record ALF2 "ALF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF2);
     annotation (preferredView = "info");
    end ALF2;

    record ALF2minus "ALF2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF2minus,
      z=-1);
     annotation (preferredView = "info");
    end ALF2minus;

    record ALF2CL "ALF2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF2CL);
     annotation (preferredView = "info");
    end ALF2CL;

    record ALF3 "ALF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF3);
     annotation (preferredView = "info");
    end ALF3;

    record ALF4minus "ALF4minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALF4minus,
      z=-1);
     annotation (preferredView = "info");
    end ALF4minus;

    record ALH "ALH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH);
     annotation (preferredView = "info");
    end ALH;

    record ALHCL "ALHCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHCL);
     annotation (preferredView = "info");
    end ALHCL;

    record ALHCL2 "ALHCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHCL2);
     annotation (preferredView = "info");
    end ALHCL2;

    record ALHF "ALHF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHF);
     annotation (preferredView = "info");
    end ALHF;

    record ALHFCL "ALHFCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHFCL);
     annotation (preferredView = "info");
    end ALHFCL;

    record ALHF2 "ALHF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALHF2);
     annotation (preferredView = "info");
    end ALHF2;

    record ALH2 "ALH2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH2);
     annotation (preferredView = "info");
    end ALH2;

    record ALH2CL "ALH2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH2CL);
     annotation (preferredView = "info");
    end ALH2CL;

    record ALH2F "ALH2F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH2F);
     annotation (preferredView = "info");
    end ALH2F;

    record ALH3 "ALH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALH3);
     annotation (preferredView = "info");
    end ALH3;

    record ALI "ALI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALI);
     annotation (preferredView = "info");
    end ALI;

    record ALI2 "ALI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALI2);
     annotation (preferredView = "info");
    end ALI2;

    record ALI3 "ALI3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALI3);
     annotation (preferredView = "info");
    end ALI3;

    record ALN "ALN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALN);
     annotation (preferredView = "info");
    end ALN;

    record ALO "ALO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALO);
     annotation (preferredView = "info");
    end ALO;

    record ALOplus "ALOplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOplus,
      z=1);
     annotation (preferredView = "info");
    end ALOplus;

    record ALOminus "ALOminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOminus,
      z=-1);
     annotation (preferredView = "info");
    end ALOminus;

    record ALOCL "ALOCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOCL);
     annotation (preferredView = "info");
    end ALOCL;

    record ALOCL2 "ALOCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOCL2);
     annotation (preferredView = "info");
    end ALOCL2;

    record ALOF "ALOF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOF);
     annotation (preferredView = "info");
    end ALOF;

    record ALOF2 "ALOF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOF2);
     annotation (preferredView = "info");
    end ALOF2;

    record ALOF2minus "ALOF2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOF2minus,
      z=-1);
     annotation (preferredView = "info");
    end ALOF2minus;

    record ALOH "ALOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOH);
     annotation (preferredView = "info");
    end ALOH;

    record ALOHCL "ALOHCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOHCL);
     annotation (preferredView = "info");
    end ALOHCL;

    record ALOHCL2 "ALOHCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOHCL2);
     annotation (preferredView = "info");
    end ALOHCL2;

    record ALOHF "ALOHF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOHF);
     annotation (preferredView = "info");
    end ALOHF;

    record ALOHF2 "ALOHF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALOHF2);
     annotation (preferredView = "info");
    end ALOHF2;

    record ALO2 "ALO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALO2);
     annotation (preferredView = "info");
    end ALO2;

    record ALO2minus "ALO2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALO2minus,
      z=-1);
     annotation (preferredView = "info");
    end ALO2minus;

    record AL_OH_2 "AL_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL_OH_2);
     annotation (preferredView = "info");
    end AL_OH_2;

    record AL_OH_2CL "AL_OH_2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL_OH_2CL);
     annotation (preferredView = "info");
    end AL_OH_2CL;

    record AL_OH_2F "AL_OH_2F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL_OH_2F);
     annotation (preferredView = "info");
    end AL_OH_2F;

    record AL_OH_3 "AL_OH_3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL_OH_3);
     annotation (preferredView = "info");
    end AL_OH_3;

    record ALS "ALS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALS);
     annotation (preferredView = "info");
    end ALS;

    record ALS2 "ALS2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ALS2);
     annotation (preferredView = "info");
    end ALS2;

    record AL2 "AL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2);
     annotation (preferredView = "info");
    end AL2;

    record AL2Br6 "AL2Br6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2Br6);
     annotation (preferredView = "info");
    end AL2Br6;

    record AL2C2 "AL2C2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2C2);
     annotation (preferredView = "info");
    end AL2C2;

    record AL2CL6 "AL2CL6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2CL6);
     annotation (preferredView = "info");
    end AL2CL6;

    record AL2F6 "AL2F6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2F6);
     annotation (preferredView = "info");
    end AL2F6;

    record AL2I6 "AL2I6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2I6);
     annotation (preferredView = "info");
    end AL2I6;

    record AL2O "AL2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2O);
     annotation (preferredView = "info");
    end AL2O;

    record AL2Oplus "AL2Oplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2Oplus,
      z=1);
     annotation (preferredView = "info");
    end AL2Oplus;

    record AL2O2 "AL2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2O2);
     annotation (preferredView = "info");
    end AL2O2;

    record AL2O2plus "AL2O2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2O2plus,
      z=1);
     annotation (preferredView = "info");
    end AL2O2plus;

    record AL2O3 "AL2O3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2O3);
     annotation (preferredView = "info");
    end AL2O3;

    record AL2S "AL2S(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2S);
     annotation (preferredView = "info");
    end AL2S;

    record AL2S2 "AL2S2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.AL2S2);
     annotation (preferredView = "info");
    end AL2S2;

    record Ar "Ar(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ar);
     annotation (preferredView = "info");
    end Ar;

    record Arplus "Arplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Arplus,
      z=1);
     annotation (preferredView = "info");
    end Arplus;

    record B "B(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B);
     annotation (preferredView = "info");
    end B;

    record Bplus "Bplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Bplus,
      z=1);
     annotation (preferredView = "info");
    end Bplus;

    record Bminus "Bminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Bminus,
      z=-1);
     annotation (preferredView = "info");
    end Bminus;

    record BBr "BBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BBr);
     annotation (preferredView = "info");
    end BBr;

    record BBr2 "BBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BBr2);
     annotation (preferredView = "info");
    end BBr2;

    record BBr3 "BBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BBr3);
     annotation (preferredView = "info");
    end BBr3;

    record BC "BC(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BC);
     annotation (preferredView = "info");
    end BC;

    record BC2 "BC2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BC2);
     annotation (preferredView = "info");
    end BC2;

    record BCL "BCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL);
     annotation (preferredView = "info");
    end BCL;

    record BCLplus "BCLplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BCLplus,
      z=1);
     annotation (preferredView = "info");
    end BCLplus;

    record BCLOH "BCLOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BCLOH);
     annotation (preferredView = "info");
    end BCLOH;

    record BCL_OH_2 "BCL_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL_OH_2);
     annotation (preferredView = "info");
    end BCL_OH_2;

    record BCL2 "BCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL2);
     annotation (preferredView = "info");
    end BCL2;

    record BCL2plus "BCL2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL2plus,
      z=1);
     annotation (preferredView = "info");
    end BCL2plus;

    record BCL2OH "BCL2OH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BCL2OH);
     annotation (preferredView = "info");
    end BCL2OH;

    record BF "BF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BF);
     annotation (preferredView = "info");
    end BF;

    record BFCL "BFCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BFCL);
     annotation (preferredView = "info");
    end BFCL;

    record BFCL2 "BFCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BFCL2);
     annotation (preferredView = "info");
    end BFCL2;

    record BFOH "BFOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BFOH);
     annotation (preferredView = "info");
    end BFOH;

    record BF_OH_2 "BF_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BF_OH_2);
     annotation (preferredView = "info");
    end BF_OH_2;

    record BF2 "BF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2);
     annotation (preferredView = "info");
    end BF2;

    record BF2plus "BF2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2plus,
      z=1);
     annotation (preferredView = "info");
    end BF2plus;

    record BF2minus "BF2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2minus,
      z=-1);
     annotation (preferredView = "info");
    end BF2minus;

    record BF2CL "BF2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2CL);
     annotation (preferredView = "info");
    end BF2CL;

    record BF2OH "BF2OH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BF2OH);
     annotation (preferredView = "info");
    end BF2OH;

    record BF3 "BF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BF3);
     annotation (preferredView = "info");
    end BF3;

    record BF4minus "BF4minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BF4minus,
      z=-1);
     annotation (preferredView = "info");
    end BF4minus;

    record BH "BH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BH);
     annotation (preferredView = "info");
    end BH;

    record BHCL "BHCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BHCL);
     annotation (preferredView = "info");
    end BHCL;

    record BHCL2 "BHCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BHCL2);
     annotation (preferredView = "info");
    end BHCL2;

    record BHF "BHF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BHF);
     annotation (preferredView = "info");
    end BHF;

    record BHFCL "BHFCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BHFCL);
     annotation (preferredView = "info");
    end BHFCL;

    record BHF2 "BHF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BHF2);
     annotation (preferredView = "info");
    end BHF2;

    record BH2 "BH2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BH2);
     annotation (preferredView = "info");
    end BH2;

    record BH2CL "BH2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BH2CL);
     annotation (preferredView = "info");
    end BH2CL;

    record BH2F "BH2F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BH2F);
     annotation (preferredView = "info");
    end BH2F;

    record BH3 "BH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BH3);
     annotation (preferredView = "info");
    end BH3;

    record BH3NH3 "BH3NH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BH3NH3);
     annotation (preferredView = "info");
    end BH3NH3;

    record BH4 "BH4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BH4);
     annotation (preferredView = "info");
    end BH4;

    record BI "BI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BI);
     annotation (preferredView = "info");
    end BI;

    record BI2 "BI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BI2);
     annotation (preferredView = "info");
    end BI2;

    record BI3 "BI3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BI3);
     annotation (preferredView = "info");
    end BI3;

    record BN "BN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BN);
     annotation (preferredView = "info");
    end BN;

    record BO "BO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BO);
     annotation (preferredView = "info");
    end BO;

    record BOminus "BOminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BOminus,
      z=-1);
     annotation (preferredView = "info");
    end BOminus;

    record BOCL "BOCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BOCL);
     annotation (preferredView = "info");
    end BOCL;

    record BOCL2 "BOCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BOCL2);
     annotation (preferredView = "info");
    end BOCL2;

    record BOF "BOF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BOF);
     annotation (preferredView = "info");
    end BOF;

    record BOF2 "BOF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BOF2);
     annotation (preferredView = "info");
    end BOF2;

    record BOH "BOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BOH);
     annotation (preferredView = "info");
    end BOH;

    record BO2 "BO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BO2);
     annotation (preferredView = "info");
    end BO2;

    record BO2minus "BO2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BO2minus,
      z=-1);
     annotation (preferredView = "info");
    end BO2minus;

    record B_OH_2 "B_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B_OH_2);
     annotation (preferredView = "info");
    end B_OH_2;

    record BS "BS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BS);
     annotation (preferredView = "info");
    end BS;

    record BS2 "BS2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BS2);
     annotation (preferredView = "info");
    end BS2;

    record B2 "B2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2);
     annotation (preferredView = "info");
    end B2;

    record B2C "B2C(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2C);
     annotation (preferredView = "info");
    end B2C;

    record B2CL4 "B2CL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2CL4);
     annotation (preferredView = "info");
    end B2CL4;

    record B2F4 "B2F4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2F4);
     annotation (preferredView = "info");
    end B2F4;

    record B2H "B2H(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H);
     annotation (preferredView = "info");
    end B2H;

    record B2H2 "B2H2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H2);
     annotation (preferredView = "info");
    end B2H2;

    record B2H3 "B2H3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H3);
     annotation (preferredView = "info");
    end B2H3;

    record B2H3_db "B2H3_db(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H3_db);
     annotation (preferredView = "info");
    end B2H3_db;

    record B2H4 "B2H4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H4);
     annotation (preferredView = "info");
    end B2H4;

    record B2H4_db "B2H4_db(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H4_db);
     annotation (preferredView = "info");
    end B2H4_db;

    record B2H5 "B2H5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H5);
     annotation (preferredView = "info");
    end B2H5;

    record B2H5_db "B2H5_db(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H5_db);
     annotation (preferredView = "info");
    end B2H5_db;

    record B2H6 "B2H6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2H6);
     annotation (preferredView = "info");
    end B2H6;

    record B2O "B2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2O);
     annotation (preferredView = "info");
    end B2O;

    record B2O2 "B2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2O2);
     annotation (preferredView = "info");
    end B2O2;

    record B2O3 "B2O3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2O3);
     annotation (preferredView = "info");
    end B2O3;

    record B2_OH_4 "B2_OH_4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2_OH_4);
     annotation (preferredView = "info");
    end B2_OH_4;

    record B2S "B2S(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2S);
     annotation (preferredView = "info");
    end B2S;

    record B2S2 "B2S2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2S2);
     annotation (preferredView = "info");
    end B2S2;

    record B2S3 "B2S3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B2S3);
     annotation (preferredView = "info");
    end B2S3;

    record B3H7_C2v "B3H7_C2v(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B3H7_C2v);
     annotation (preferredView = "info");
    end B3H7_C2v;

    record B3H7_Cs "B3H7_Cs(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B3H7_Cs);
     annotation (preferredView = "info");
    end B3H7_Cs;

    record B3H9 "B3H9(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B3H9);
     annotation (preferredView = "info");
    end B3H9;

    record B3N3H6 "B3N3H6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B3N3H6);
     annotation (preferredView = "info");
    end B3N3H6;

    record B3O3CL3 "B3O3CL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B3O3CL3);
     annotation (preferredView = "info");
    end B3O3CL3;

    record B3O3FCL2 "B3O3FCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B3O3FCL2);
     annotation (preferredView = "info");
    end B3O3FCL2;

    record B3O3F2CL "B3O3F2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B3O3F2CL);
     annotation (preferredView = "info");
    end B3O3F2CL;

    record B3O3F3 "B3O3F3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B3O3F3);
     annotation (preferredView = "info");
    end B3O3F3;

    record B4H4 "B4H4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B4H4);
     annotation (preferredView = "info");
    end B4H4;

    record B4H10 "B4H10(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B4H10);
     annotation (preferredView = "info");
    end B4H10;

    record B4H12 "B4H12(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B4H12);
     annotation (preferredView = "info");
    end B4H12;

    record B5H9 "B5H9(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.B5H9);
     annotation (preferredView = "info");
    end B5H9;

    record Ba "Ba(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ba);
     annotation (preferredView = "info");
    end Ba;

    record Baplus "Baplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Baplus,
      z=1);
     annotation (preferredView = "info");
    end Baplus;

    record BaBr "BaBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaBr);
     annotation (preferredView = "info");
    end BaBr;

    record BaBr2 "BaBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaBr2);
     annotation (preferredView = "info");
    end BaBr2;

    record BaCL "BaCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaCL);
     annotation (preferredView = "info");
    end BaCL;

    record BaCLplus "BaCLplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaCLplus,
      z=1);
     annotation (preferredView = "info");
    end BaCLplus;

    record BaCL2 "BaCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaCL2);
     annotation (preferredView = "info");
    end BaCL2;

    record BaF "BaF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaF);
     annotation (preferredView = "info");
    end BaF;

    record BaFplus "BaFplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaFplus,
      z=1);
     annotation (preferredView = "info");
    end BaFplus;

    record BaF2 "BaF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaF2);
     annotation (preferredView = "info");
    end BaF2;

    record BaH "BaH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaH);
     annotation (preferredView = "info");
    end BaH;

    record BaI "BaI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaI);
     annotation (preferredView = "info");
    end BaI;

    record BaI2 "BaI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaI2);
     annotation (preferredView = "info");
    end BaI2;

    record BaO "BaO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaO);
     annotation (preferredView = "info");
    end BaO;

    record BaOplus "BaOplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaOplus,
      z=1);
     annotation (preferredView = "info");
    end BaOplus;

    record BaOH "BaOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaOH);
     annotation (preferredView = "info");
    end BaOH;

    record BaOHplus "BaOHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaOHplus,
      z=1);
     annotation (preferredView = "info");
    end BaOHplus;

    record Ba_OH_2 "Ba_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ba_OH_2);
     annotation (preferredView = "info");
    end Ba_OH_2;

    record BaS "BaS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BaS);
     annotation (preferredView = "info");
    end BaS;

    record Ba2 "Ba2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ba2);
     annotation (preferredView = "info");
    end Ba2;

    record Be "Be(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be);
     annotation (preferredView = "info");
    end Be;

    record Beplus "Beplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Beplus,
      z=1);
     annotation (preferredView = "info");
    end Beplus;

    record Beplusplus "Beplusplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Beplusplus,
      z=1);
     annotation (preferredView = "info");
    end Beplusplus;

    record BeBr "BeBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeBr);
     annotation (preferredView = "info");
    end BeBr;

    record BeBr2 "BeBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeBr2);
     annotation (preferredView = "info");
    end BeBr2;

    record BeCL "BeCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeCL);
     annotation (preferredView = "info");
    end BeCL;

    record BeCL2 "BeCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeCL2);
     annotation (preferredView = "info");
    end BeCL2;

    record BeF "BeF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeF);
     annotation (preferredView = "info");
    end BeF;

    record BeF2 "BeF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeF2);
     annotation (preferredView = "info");
    end BeF2;

    record BeH "BeH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeH);
     annotation (preferredView = "info");
    end BeH;

    record BeHplus "BeHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeHplus,
      z=1);
     annotation (preferredView = "info");
    end BeHplus;

    record BeH2 "BeH2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeH2);
     annotation (preferredView = "info");
    end BeH2;

    record BeI "BeI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeI);
     annotation (preferredView = "info");
    end BeI;

    record BeI2 "BeI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeI2);
     annotation (preferredView = "info");
    end BeI2;

    record BeN "BeN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeN);
     annotation (preferredView = "info");
    end BeN;

    record BeO "BeO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeO);
     annotation (preferredView = "info");
    end BeO;

    record BeOH "BeOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeOH);
     annotation (preferredView = "info");
    end BeOH;

    record BeOHplus "BeOHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeOHplus,
      z=1);
     annotation (preferredView = "info");
    end BeOHplus;

    record Be_OH_2 "Be_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be_OH_2);
     annotation (preferredView = "info");
    end Be_OH_2;

    record BeS "BeS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BeS);
     annotation (preferredView = "info");
    end BeS;

    record Be2 "Be2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2);
     annotation (preferredView = "info");
    end Be2;

    record Be2CL4 "Be2CL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2CL4);
     annotation (preferredView = "info");
    end Be2CL4;

    record Be2F4 "Be2F4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2F4);
     annotation (preferredView = "info");
    end Be2F4;

    record Be2O "Be2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2O);
     annotation (preferredView = "info");
    end Be2O;

    record Be2OF2 "Be2OF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2OF2);
     annotation (preferredView = "info");
    end Be2OF2;

    record Be2O2 "Be2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be2O2);
     annotation (preferredView = "info");
    end Be2O2;

    record Be3O3 "Be3O3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be3O3);
     annotation (preferredView = "info");
    end Be3O3;

    record Be4O4 "Be4O4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Be4O4);
     annotation (preferredView = "info");
    end Be4O4;

    record Br "Br(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Br);
     annotation (preferredView = "info");
    end Br;

    record Brplus "Brplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Brplus,
      z=1);
     annotation (preferredView = "info");
    end Brplus;

    record Brminus "Brminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Brminus,
      z=-1);
     annotation (preferredView = "info");
    end Brminus;

    record BrCL "BrCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BrCL);
     annotation (preferredView = "info");
    end BrCL;

    record BrF "BrF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BrF);
     annotation (preferredView = "info");
    end BrF;

    record BrF3 "BrF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BrF3);
     annotation (preferredView = "info");
    end BrF3;

    record BrF5 "BrF5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BrF5);
     annotation (preferredView = "info");
    end BrF5;

    record BrO "BrO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BrO);
     annotation (preferredView = "info");
    end BrO;

    record OBrO "OBrO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.OBrO);
     annotation (preferredView = "info");
    end OBrO;

    record BrOO "BrOO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BrOO);
     annotation (preferredView = "info");
    end BrOO;

    record BrO3 "BrO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BrO3);
     annotation (preferredView = "info");
    end BrO3;

    record Br2 "Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Br2);
     annotation (preferredView = "info");
    end Br2;

    record BrBrO "BrBrO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BrBrO);
     annotation (preferredView = "info");
    end BrBrO;

    record BrOBr "BrOBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.BrOBr);
     annotation (preferredView = "info");
    end BrOBr;

    record C "C(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C);
     annotation (preferredView = "info");
    end C;

    record Cplus "Cplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cplus,
      z=1);
     annotation (preferredView = "info");
    end Cplus;

    record Cminus "Cminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cminus,
      z=-1);
     annotation (preferredView = "info");
    end Cminus;

    record CBr "CBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CBr);
     annotation (preferredView = "info");
    end CBr;

    record CBr2 "CBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CBr2);
     annotation (preferredView = "info");
    end CBr2;

    record CBr3 "CBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CBr3);
     annotation (preferredView = "info");
    end CBr3;

    record CBr4 "CBr4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CBr4);
     annotation (preferredView = "info");
    end CBr4;

    record CCL "CCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL);
     annotation (preferredView = "info");
    end CCL;

    record CCL2 "CCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL2);
     annotation (preferredView = "info");
    end CCL2;

    record CCL2Br2 "CCL2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL2Br2);
     annotation (preferredView = "info");
    end CCL2Br2;

    record CCL3 "CCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL3);
     annotation (preferredView = "info");
    end CCL3;

    record CCL3Br "CCL3Br(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL3Br);
     annotation (preferredView = "info");
    end CCL3Br;

    record CCL4 "CCL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CCL4);
     annotation (preferredView = "info");
    end CCL4;

    record CF "CF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF);
     annotation (preferredView = "info");
    end CF;

    record CFplus "CFplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CFplus,
      z=1);
     annotation (preferredView = "info");
    end CFplus;

    record CFBr3 "CFBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CFBr3);
     annotation (preferredView = "info");
    end CFBr3;

    record CFCL "CFCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCL);
     annotation (preferredView = "info");
    end CFCL;

    record CFCLBr2 "CFCLBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCLBr2);
     annotation (preferredView = "info");
    end CFCLBr2;

    record CFCL2 "CFCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCL2);
     annotation (preferredView = "info");
    end CFCL2;

    record CFCL2Br "CFCL2Br(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCL2Br);
     annotation (preferredView = "info");
    end CFCL2Br;

    record CFCL3 "CFCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CFCL3);
     annotation (preferredView = "info");
    end CFCL3;

    record CF2 "CF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2);
     annotation (preferredView = "info");
    end CF2;

    record CF2plus "CF2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2plus,
      z=1);
     annotation (preferredView = "info");
    end CF2plus;

    record CF2Br2 "CF2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2Br2);
     annotation (preferredView = "info");
    end CF2Br2;

    record CF2CL "CF2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2CL);
     annotation (preferredView = "info");
    end CF2CL;

    record CF2CLBr "CF2CLBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2CLBr);
     annotation (preferredView = "info");
    end CF2CLBr;

    record CF2CL2 "CF2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF2CL2);
     annotation (preferredView = "info");
    end CF2CL2;

    record CF3 "CF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF3);
     annotation (preferredView = "info");
    end CF3;

    record CF3plus "CF3plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF3plus,
      z=1);
     annotation (preferredView = "info");
    end CF3plus;

    record CF3Br "CF3Br(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF3Br);
     annotation (preferredView = "info");
    end CF3Br;

    record CF3CL "CF3CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF3CL);
     annotation (preferredView = "info");
    end CF3CL;

    record CF4 "CF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CF4);
     annotation (preferredView = "info");
    end CF4;

    record CHplus "CHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHplus,
      z=1);
     annotation (preferredView = "info");
    end CHplus;

    record CHBr3 "CHBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHBr3);
     annotation (preferredView = "info");
    end CHBr3;

    record CHCL "CHCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCL);
     annotation (preferredView = "info");
    end CHCL;

    record CHCLBr2 "CHCLBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCLBr2);
     annotation (preferredView = "info");
    end CHCLBr2;

    record CHCL2 "CHCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCL2);
     annotation (preferredView = "info");
    end CHCL2;

    record CHCL2Br "CHCL2Br(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCL2Br);
     annotation (preferredView = "info");
    end CHCL2Br;

    record CHCL3 "CHCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHCL3);
     annotation (preferredView = "info");
    end CHCL3;

    record CHF "CHF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF);
     annotation (preferredView = "info");
    end CHF;

    record CHFBr2 "CHFBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHFBr2);
     annotation (preferredView = "info");
    end CHFBr2;

    record CHFCL "CHFCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHFCL);
     annotation (preferredView = "info");
    end CHFCL;

    record CHFCLBr "CHFCLBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHFCLBr);
     annotation (preferredView = "info");
    end CHFCLBr;

    record CHFCL2 "CHFCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHFCL2);
     annotation (preferredView = "info");
    end CHFCL2;

    record CHF2 "CHF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF2);
     annotation (preferredView = "info");
    end CHF2;

    record CHF2Br "CHF2Br(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF2Br);
     annotation (preferredView = "info");
    end CHF2Br;

    record CHF2CL "CHF2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF2CL);
     annotation (preferredView = "info");
    end CHF2CL;

    record CHF3 "CHF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHF3);
     annotation (preferredView = "info");
    end CHF3;

    record CHI3 "CHI3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CHI3);
     annotation (preferredView = "info");
    end CHI3;

    record CH2 "CH2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2);
     annotation (preferredView = "info");
    end CH2;

    record CH2Br2 "CH2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2Br2);
     annotation (preferredView = "info");
    end CH2Br2;

    record CH2CL "CH2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CL);
     annotation (preferredView = "info");
    end CH2CL;

    record CH2CLBr "CH2CLBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CLBr);
     annotation (preferredView = "info");
    end CH2CLBr;

    record CH2CL2 "CH2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CL2);
     annotation (preferredView = "info");
    end CH2CL2;

    record CH2F "CH2F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2F);
     annotation (preferredView = "info");
    end CH2F;

    record CH2FBr "CH2FBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2FBr);
     annotation (preferredView = "info");
    end CH2FBr;

    record CH2FCL "CH2FCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2FCL);
     annotation (preferredView = "info");
    end CH2FCL;

    record CH2F2 "CH2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2F2);
     annotation (preferredView = "info");
    end CH2F2;

    record CH2I2 "CH2I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2I2);
     annotation (preferredView = "info");
    end CH2I2;

    record CH3 "CH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3);
     annotation (preferredView = "info");
    end CH3;

    record CH3Br "CH3Br(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3Br);
     annotation (preferredView = "info");
    end CH3Br;

    record CH3CL "CH3CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3CL);
     annotation (preferredView = "info");
    end CH3CL;

    record CH3F "CH3F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3F);
     annotation (preferredView = "info");
    end CH3F;

    record CH3I "CH3I(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3I);
     annotation (preferredView = "info");
    end CH3I;

    record CH2OH "CH2OH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2OH);
     annotation (preferredView = "info");
    end CH2OH;

    record CH2OHplus "CH2OHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2OHplus,
      z=1);
     annotation (preferredView = "info");
    end CH2OHplus;

    record CH3O "CH3O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3O);
     annotation (preferredView = "info");
    end CH3O;

    record CH4 "CH4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH4);
     annotation (preferredView = "info");
    end CH4;

    record CH3OH "CH3OH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3OH);
     annotation (preferredView = "info");
    end CH3OH;

    record CH3OOH "CH3OOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3OOH);
     annotation (preferredView = "info");
    end CH3OOH;

    record CI "CI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CI);
     annotation (preferredView = "info");
    end CI;

    record CI2 "CI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CI2);
     annotation (preferredView = "info");
    end CI2;

    record CI3 "CI3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CI3);
     annotation (preferredView = "info");
    end CI3;

    record CI4 "CI4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CI4);
     annotation (preferredView = "info");
    end CI4;

    record CN "CN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CN);
     annotation (preferredView = "info");
    end CN;

    record CNplus "CNplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CNplus,
      z=1);
     annotation (preferredView = "info");
    end CNplus;

    record CNminus "CNminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CNminus,
      z=-1);
     annotation (preferredView = "info");
    end CNminus;

    record CNN "CNN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CNN);
     annotation (preferredView = "info");
    end CNN;

    record CO "CO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CO);
     annotation (preferredView = "info");
    end CO;

    record COplus "COplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.COplus,
      z=1);
     annotation (preferredView = "info");
    end COplus;

    record COCL "COCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.COCL);
     annotation (preferredView = "info");
    end COCL;

    record COCL2 "COCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.COCL2);
     annotation (preferredView = "info");
    end COCL2;

    record COFCL "COFCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.COFCL);
     annotation (preferredView = "info");
    end COFCL;

    record COF2 "COF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.COF2);
     annotation (preferredView = "info");
    end COF2;

    record COHCL "COHCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.COHCL);
     annotation (preferredView = "info");
    end COHCL;

    record COHF "COHF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.COHF);
     annotation (preferredView = "info");
    end COHF;

    record COS "COS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.COS);
     annotation (preferredView = "info");
    end COS;

    record CO2 "CO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CO2);
     annotation (preferredView = "info");
    end CO2;

    record CO2plus "CO2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CO2plus,
      z=1);
     annotation (preferredView = "info");
    end CO2plus;

    record COOH "COOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.COOH);
     annotation (preferredView = "info");
    end COOH;

    record CP "CP(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CP);
     annotation (preferredView = "info");
    end CP;

    record CS "CS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CS);
     annotation (preferredView = "info");
    end CS;

    record CS2 "CS2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CS2);
     annotation (preferredView = "info");
    end CS2;

    record C2 "C2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2);
     annotation (preferredView = "info");
    end C2;

    record C2plus "C2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2plus,
      z=1);
     annotation (preferredView = "info");
    end C2plus;

    record C2minus "C2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2minus,
      z=-1);
     annotation (preferredView = "info");
    end C2minus;

    record C2CL "C2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL);
     annotation (preferredView = "info");
    end C2CL;

    record C2CL2 "C2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL2);
     annotation (preferredView = "info");
    end C2CL2;

    record C2CL3 "C2CL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL3);
     annotation (preferredView = "info");
    end C2CL3;

    record C2CL4 "C2CL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL4);
     annotation (preferredView = "info");
    end C2CL4;

    record C2CL6 "C2CL6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2CL6);
     annotation (preferredView = "info");
    end C2CL6;

    record C2F "C2F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F);
     annotation (preferredView = "info");
    end C2F;

    record C2FCL "C2FCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2FCL);
     annotation (preferredView = "info");
    end C2FCL;

    record C2FCL3 "C2FCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2FCL3);
     annotation (preferredView = "info");
    end C2FCL3;

    record C2F2 "C2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F2);
     annotation (preferredView = "info");
    end C2F2;

    record C2F2CL2 "C2F2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F2CL2);
     annotation (preferredView = "info");
    end C2F2CL2;

    record C2F3 "C2F3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F3);
     annotation (preferredView = "info");
    end C2F3;

    record C2F3CL "C2F3CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F3CL);
     annotation (preferredView = "info");
    end C2F3CL;

    record C2F4 "C2F4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F4);
     annotation (preferredView = "info");
    end C2F4;

    record C2F6 "C2F6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2F6);
     annotation (preferredView = "info");
    end C2F6;

    record C2H "C2H(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H);
     annotation (preferredView = "info");
    end C2H;

    record C2HCL "C2HCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HCL);
     annotation (preferredView = "info");
    end C2HCL;

    record C2HCL3 "C2HCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HCL3);
     annotation (preferredView = "info");
    end C2HCL3;

    record C2HF "C2HF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HF);
     annotation (preferredView = "info");
    end C2HF;

    record C2HFCL2 "C2HFCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HFCL2);
     annotation (preferredView = "info");
    end C2HFCL2;

    record C2HF2CL "C2HF2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HF2CL);
     annotation (preferredView = "info");
    end C2HF2CL;

    record C2HF3 "C2HF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2HF3);
     annotation (preferredView = "info");
    end C2HF3;

    record C2H2_vinylidene "C2H2_vinylidene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H2_vinylidene);
     annotation (preferredView = "info");
    end C2H2_vinylidene;

    record C2H2CL2 "C2H2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H2CL2);
     annotation (preferredView = "info");
    end C2H2CL2;

    record C2H2FCL "C2H2FCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H2FCL);
     annotation (preferredView = "info");
    end C2H2FCL;

    record C2H2F2 "C2H2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H2F2);
     annotation (preferredView = "info");
    end C2H2F2;

    record CH2CO_ketene "CH2CO_ketene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CO_ketene);
     annotation (preferredView = "info");
    end CH2CO_ketene;

    record O_CH_2O "O_CH_2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.O_CH_2O);
     annotation (preferredView = "info");
    end O_CH_2O;

    record HO_CO_2OH "HO_CO_2OH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HO_CO_2OH);
     annotation (preferredView = "info");
    end HO_CO_2OH;

    record C2H3_vinyl "C2H3_vinyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H3_vinyl);
     annotation (preferredView = "info");
    end C2H3_vinyl;

    record CH2BrminusCOOH "CH2BrminusCOOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2BrminusCOOH);
     annotation (preferredView = "info");
    end CH2BrminusCOOH;

    record C2H3CL "C2H3CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H3CL);
     annotation (preferredView = "info");
    end C2H3CL;

    record CH2CLminusCOOH "CH2CLminusCOOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH2CLminusCOOH);
     annotation (preferredView = "info");
    end CH2CLminusCOOH;

    record C2H3F "C2H3F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H3F);
     annotation (preferredView = "info");
    end C2H3F;

    record CH3CN "CH3CN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3CN);
     annotation (preferredView = "info");
    end CH3CN;

    record CH3CO_acetyl "CH3CO_acetyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3CO_acetyl);
     annotation (preferredView = "info");
    end CH3CO_acetyl;

    record C2H4 "C2H4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H4);
     annotation (preferredView = "info");
    end C2H4;

    record C2H4O_ethylen_o "C2H4O_ethylen_o(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H4O_ethylen_o);
     annotation (preferredView = "info");
    end C2H4O_ethylen_o;

    record CH3CHO_ethanal "CH3CHO_ethanal(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3CHO_ethanal);
     annotation (preferredView = "info");
    end CH3CHO_ethanal;

    record CH3COOH "CH3COOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3COOH);
     annotation (preferredView = "info");
    end CH3COOH;

    record OHCH2COOH "OHCH2COOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.OHCH2COOH);
     annotation (preferredView = "info");
    end OHCH2COOH;

    record C2H5 "C2H5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H5);
     annotation (preferredView = "info");
    end C2H5;

    record C2H5Br "C2H5Br(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H5Br);
     annotation (preferredView = "info");
    end C2H5Br;

    record C2H6 "C2H6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H6);
     annotation (preferredView = "info");
    end C2H6;

    record CH3N2CH3 "CH3N2CH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3N2CH3);
     annotation (preferredView = "info");
    end CH3N2CH3;

    record C2H5OH "C2H5OH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2H5OH);
     annotation (preferredView = "info");
    end C2H5OH;

    record CH3OCH3 "CH3OCH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3OCH3);
     annotation (preferredView = "info");
    end CH3OCH3;

    record CH3O2CH3 "CH3O2CH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3O2CH3);
     annotation (preferredView = "info");
    end CH3O2CH3;

    record CCN "CCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CCN);
     annotation (preferredView = "info");
    end CCN;

    record CNC "CNC(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CNC);
     annotation (preferredView = "info");
    end CNC;

    record OCCN "OCCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.OCCN);
     annotation (preferredView = "info");
    end OCCN;

    record C2N2 "C2N2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2N2);
     annotation (preferredView = "info");
    end C2N2;

    record C2O "C2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C2O);
     annotation (preferredView = "info");
    end C2O;

    record C3 "C3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3);
     annotation (preferredView = "info");
    end C3;

    record C3H3_1_propynl "C3H3_1_propynl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H3_1_propynl);
     annotation (preferredView = "info");
    end C3H3_1_propynl;

    record C3H3_2_propynl "C3H3_2_propynl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H3_2_propynl);
     annotation (preferredView = "info");
    end C3H3_2_propynl;

    record C3H4_allene "C3H4_allene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H4_allene);
     annotation (preferredView = "info");
    end C3H4_allene;

    record C3H4_propyne "C3H4_propyne(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H4_propyne);
     annotation (preferredView = "info");
    end C3H4_propyne;

    record C3H4_cyclo "C3H4_cyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H4_cyclo);
     annotation (preferredView = "info");
    end C3H4_cyclo;

    record C3H5_allyl "C3H5_allyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H5_allyl);
     annotation (preferredView = "info");
    end C3H5_allyl;

    record C3H6_propylene "C3H6_propylene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6_propylene);
     annotation (preferredView = "info");
    end C3H6_propylene;

    record C3H6_cyclo "C3H6_cyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6_cyclo);
     annotation (preferredView = "info");
    end C3H6_cyclo;

    record C3H6O_propylox "C3H6O_propylox(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6O_propylox);
     annotation (preferredView = "info");
    end C3H6O_propylox;

    record C3H6O_acetone "C3H6O_acetone(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6O_acetone);
     annotation (preferredView = "info");
    end C3H6O_acetone;

    record C3H6O_propanal "C3H6O_propanal(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H6O_propanal);
     annotation (preferredView = "info");
    end C3H6O_propanal;

    record C3H7_n_propyl "C3H7_n_propyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H7_n_propyl);
     annotation (preferredView = "info");
    end C3H7_n_propyl;

    record C3H7_i_propyl "C3H7_i_propyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H7_i_propyl);
     annotation (preferredView = "info");
    end C3H7_i_propyl;

    record C3H8 "C3H8(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H8);
     annotation (preferredView = "info");
    end C3H8;

    record C3H8O_1propanol "C3H8O_1propanol(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H8O_1propanol);
     annotation (preferredView = "info");
    end C3H8O_1propanol;

    record C3H8O_2propanol "C3H8O_2propanol(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3H8O_2propanol);
     annotation (preferredView = "info");
    end C3H8O_2propanol;

    record CNCOCN "CNCOCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CNCOCN);
     annotation (preferredView = "info");
    end CNCOCN;

    record C3O2 "C3O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C3O2);
     annotation (preferredView = "info");
    end C3O2;

    record C4 "C4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4);
     annotation (preferredView = "info");
    end C4;

    record C4H2_butadiyne "C4H2_butadiyne(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H2_butadiyne);
     annotation (preferredView = "info");
    end C4H2_butadiyne;

    record C4H4_1_3minuscyclo "C4H4_1_3minuscyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H4_1_3minuscyclo);
     annotation (preferredView = "info");
    end C4H4_1_3minuscyclo;

    record C4H6_butadiene "C4H6_butadiene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H6_butadiene);
     annotation (preferredView = "info");
    end C4H6_butadiene;

    record C4H6_1butyne "C4H6_1butyne(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H6_1butyne);
     annotation (preferredView = "info");
    end C4H6_1butyne;

    record C4H6_2butyne "C4H6_2butyne(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H6_2butyne);
     annotation (preferredView = "info");
    end C4H6_2butyne;

    record C4H6_cyclo "C4H6_cyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H6_cyclo);
     annotation (preferredView = "info");
    end C4H6_cyclo;

    record C4H8_1_butene "C4H8_1_butene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H8_1_butene);
     annotation (preferredView = "info");
    end C4H8_1_butene;

    record C4H8_cis2_buten "C4H8_cis2_buten(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H8_cis2_buten);
     annotation (preferredView = "info");
    end C4H8_cis2_buten;

    record C4H8_isobutene "C4H8_isobutene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H8_isobutene);
     annotation (preferredView = "info");
    end C4H8_isobutene;

    record C4H8_cyclo "C4H8_cyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H8_cyclo);
     annotation (preferredView = "info");
    end C4H8_cyclo;

    record C4H9_n_butyl "C4H9_n_butyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H9_n_butyl);
     annotation (preferredView = "info");
    end C4H9_n_butyl;

    record C4H9_i_butyl "C4H9_i_butyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H9_i_butyl);
     annotation (preferredView = "info");
    end C4H9_i_butyl;

    record C4H9_s_butyl "C4H9_s_butyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H9_s_butyl);
     annotation (preferredView = "info");
    end C4H9_s_butyl;

    record C4H9_t_butyl "C4H9_t_butyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H9_t_butyl);
     annotation (preferredView = "info");
    end C4H9_t_butyl;

    record C4H10_n_butane "C4H10_n_butane(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H10_n_butane);
     annotation (preferredView = "info");
    end C4H10_n_butane;

    record C4H10_isobutane "C4H10_isobutane(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4H10_isobutane);
     annotation (preferredView = "info");
    end C4H10_isobutane;

    record C4N2 "C4N2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C4N2);
     annotation (preferredView = "info");
    end C4N2;

    record C5 "C5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C5);
     annotation (preferredView = "info");
    end C5;

    record C5H6_1_3cyclo "C5H6_1_3cyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H6_1_3cyclo);
     annotation (preferredView = "info");
    end C5H6_1_3cyclo;

    record C5H8_cyclo "C5H8_cyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H8_cyclo);
     annotation (preferredView = "info");
    end C5H8_cyclo;

    record C5H10_1_pentene "C5H10_1_pentene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H10_1_pentene);
     annotation (preferredView = "info");
    end C5H10_1_pentene;

    record C5H10_cyclo "C5H10_cyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H10_cyclo);
     annotation (preferredView = "info");
    end C5H10_cyclo;

    record C5H11_pentyl "C5H11_pentyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H11_pentyl);
     annotation (preferredView = "info");
    end C5H11_pentyl;

    record C5H11_t_pentyl "C5H11_t_pentyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H11_t_pentyl);
     annotation (preferredView = "info");
    end C5H11_t_pentyl;

    record C5H12_n_pentane "C5H12_n_pentane(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H12_n_pentane);
     annotation (preferredView = "info");
    end C5H12_n_pentane;

    record C5H12_i_pentane "C5H12_i_pentane(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C5H12_i_pentane);
     annotation (preferredView = "info");
    end C5H12_i_pentane;

    record CH3C_CH3_2CH3 "CH3C_CH3_2CH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CH3C_CH3_2CH3);
     annotation (preferredView = "info");
    end CH3C_CH3_2CH3;

    record C6D5_phenyl "C6D5_phenyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6D5_phenyl);
     annotation (preferredView = "info");
    end C6D5_phenyl;

    record C6D6 "C6D6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6D6);
     annotation (preferredView = "info");
    end C6D6;

    record C6H2 "C6H2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H2);
     annotation (preferredView = "info");
    end C6H2;

    record C6H5_phenyl "C6H5_phenyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H5_phenyl);
     annotation (preferredView = "info");
    end C6H5_phenyl;

    record C6H5O_phenoxy "C6H5O_phenoxy(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H5O_phenoxy);
     annotation (preferredView = "info");
    end C6H5O_phenoxy;

    record C6H6 "C6H6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H6);
     annotation (preferredView = "info");
    end C6H6;

    record C6H5OH_phenol "C6H5OH_phenol(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H5OH_phenol);
     annotation (preferredView = "info");
    end C6H5OH_phenol;

    record C6H10_cyclo "C6H10_cyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H10_cyclo);
     annotation (preferredView = "info");
    end C6H10_cyclo;

    record C6H12_1_hexene "C6H12_1_hexene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H12_1_hexene);
     annotation (preferredView = "info");
    end C6H12_1_hexene;

    record C6H12_cyclo "C6H12_cyclo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H12_cyclo);
     annotation (preferredView = "info");
    end C6H12_cyclo;

    record C6H13_n_hexyl "C6H13_n_hexyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H13_n_hexyl);
     annotation (preferredView = "info");
    end C6H13_n_hexyl;

    record C6H14_n_hexane "C6H14_n_hexane(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C6H14_n_hexane);
     annotation (preferredView = "info");
    end C6H14_n_hexane;

    record C7H7_benzyl "C7H7_benzyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H7_benzyl);
     annotation (preferredView = "info");
    end C7H7_benzyl;

    record C7H8 "C7H8(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H8);
     annotation (preferredView = "info");
    end C7H8;

    record C7H8O_cresol_mx "C7H8O_cresol_mx(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H8O_cresol_mx);
     annotation (preferredView = "info");
    end C7H8O_cresol_mx;

    record C7H14_1_heptene "C7H14_1_heptene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H14_1_heptene);
     annotation (preferredView = "info");
    end C7H14_1_heptene;

    record C7H15_n_heptyl "C7H15_n_heptyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H15_n_heptyl);
     annotation (preferredView = "info");
    end C7H15_n_heptyl;

    record C7H16_n_heptane "C7H16_n_heptane(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H16_n_heptane);
     annotation (preferredView = "info");
    end C7H16_n_heptane;

    record C7H16_2_methylh "C7H16_2_methylh(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C7H16_2_methylh);
     annotation (preferredView = "info");
    end C7H16_2_methylh;

    record C8H8_styrene "C8H8_styrene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H8_styrene);
     annotation (preferredView = "info");
    end C8H8_styrene;

    record C8H10_ethylbenz "C8H10_ethylbenz(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H10_ethylbenz);
     annotation (preferredView = "info");
    end C8H10_ethylbenz;

    record C8H16_1_octene "C8H16_1_octene(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H16_1_octene);
     annotation (preferredView = "info");
    end C8H16_1_octene;

    record C8H17_n_octyl "C8H17_n_octyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H17_n_octyl);
     annotation (preferredView = "info");
    end C8H17_n_octyl;

    record C8H18_n_octane "C8H18_n_octane(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H18_n_octane);
     annotation (preferredView = "info");
    end C8H18_n_octane;

    record C8H18_isooctane "C8H18_isooctane(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C8H18_isooctane);
     annotation (preferredView = "info");
    end C8H18_isooctane;

    record C9H19_n_nonyl "C9H19_n_nonyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C9H19_n_nonyl);
     annotation (preferredView = "info");
    end C9H19_n_nonyl;

    record C10H8_naphthale "C10H8_naphthale(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C10H8_naphthale);
     annotation (preferredView = "info");
    end C10H8_naphthale;

    record C10H21_n_decyl "C10H21_n_decyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C10H21_n_decyl);
     annotation (preferredView = "info");
    end C10H21_n_decyl;

    record C12H9_o_bipheny "C12H9_o_bipheny(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C12H9_o_bipheny);
     annotation (preferredView = "info");
    end C12H9_o_bipheny;

    record C12H10_biphenyl "C12H10_biphenyl(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.C12H10_biphenyl);
     annotation (preferredView = "info");
    end C12H10_biphenyl;

    record Ca "Ca(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ca);
     annotation (preferredView = "info");
    end Ca;

    record Caplus "Caplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Caplus,
      z=1);
     annotation (preferredView = "info");
    end Caplus;

    record CaBr "CaBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaBr);
     annotation (preferredView = "info");
    end CaBr;

    record CaBr2 "CaBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaBr2);
     annotation (preferredView = "info");
    end CaBr2;

    record CaCL "CaCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaCL);
     annotation (preferredView = "info");
    end CaCL;

    record CaCLplus "CaCLplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaCLplus,
      z=1);
     annotation (preferredView = "info");
    end CaCLplus;

    record CaCL2 "CaCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaCL2);
     annotation (preferredView = "info");
    end CaCL2;

    record CaF "CaF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaF);
     annotation (preferredView = "info");
    end CaF;

    record CaFplus "CaFplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaFplus,
      z=1);
     annotation (preferredView = "info");
    end CaFplus;

    record CaF2 "CaF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaF2);
     annotation (preferredView = "info");
    end CaF2;

    record CaH "CaH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaH);
     annotation (preferredView = "info");
    end CaH;

    record CaI "CaI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaI);
     annotation (preferredView = "info");
    end CaI;

    record CaI2 "CaI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaI2);
     annotation (preferredView = "info");
    end CaI2;

    record CaO "CaO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaO);
     annotation (preferredView = "info");
    end CaO;

    record CaOplus "CaOplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaOplus,
      z=1);
     annotation (preferredView = "info");
    end CaOplus;

    record CaOH "CaOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaOH);
     annotation (preferredView = "info");
    end CaOH;

    record CaOHplus "CaOHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaOHplus,
      z=1);
     annotation (preferredView = "info");
    end CaOHplus;

    record Ca_OH_2 "Ca_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ca_OH_2);
     annotation (preferredView = "info");
    end Ca_OH_2;

    record CaS "CaS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CaS);
     annotation (preferredView = "info");
    end CaS;

    record Ca2 "Ca2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ca2);
     annotation (preferredView = "info");
    end Ca2;

    record Cd "Cd(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cd);
     annotation (preferredView = "info");
    end Cd;

    record Cdplus "Cdplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cdplus,
      z=1);
     annotation (preferredView = "info");
    end Cdplus;

    record CL "CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CL);
     annotation (preferredView = "info");
    end CL;

    record CLplus "CLplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CLplus,
      z=1);
     annotation (preferredView = "info");
    end CLplus;

    record CLminus "CLminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CLminus,
      z=-1);
     annotation (preferredView = "info");
    end CLminus;

    record CLCN "CLCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CLCN);
     annotation (preferredView = "info");
    end CLCN;

    record CLF "CLF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CLF);
     annotation (preferredView = "info");
    end CLF;

    record CLF3 "CLF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CLF3);
     annotation (preferredView = "info");
    end CLF3;

    record CLF5 "CLF5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CLF5);
     annotation (preferredView = "info");
    end CLF5;

    record CLO "CLO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CLO);
     annotation (preferredView = "info");
    end CLO;

    record CLO2 "CLO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CLO2);
     annotation (preferredView = "info");
    end CLO2;

    record CL2 "CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CL2);
     annotation (preferredView = "info");
    end CL2;

    record CL2O "CL2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CL2O);
     annotation (preferredView = "info");
    end CL2O;

    record Co "Co(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Co);
     annotation (preferredView = "info");
    end Co;

    record Coplus "Coplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Coplus,
      z=1);
     annotation (preferredView = "info");
    end Coplus;

    record Cominus "Cominus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cominus,
      z=-1);
     annotation (preferredView = "info");
    end Cominus;

    record Cr "Cr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cr);
     annotation (preferredView = "info");
    end Cr;

    record Crplus "Crplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Crplus,
      z=1);
     annotation (preferredView = "info");
    end Crplus;

    record Crminus "Crminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Crminus,
      z=-1);
     annotation (preferredView = "info");
    end Crminus;

    record CrN "CrN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CrN);
     annotation (preferredView = "info");
    end CrN;

    record CrO "CrO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CrO);
     annotation (preferredView = "info");
    end CrO;

    record CrO2 "CrO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CrO2);
     annotation (preferredView = "info");
    end CrO2;

    record CrO3 "CrO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CrO3);
     annotation (preferredView = "info");
    end CrO3;

    record CrO3minus "CrO3minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CrO3minus,
      z=-1);
     annotation (preferredView = "info");
    end CrO3minus;

    record Cs "Cs(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs);
     annotation (preferredView = "info");
    end Cs;

    record Csplus "Csplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Csplus,
      z=1);
     annotation (preferredView = "info");
    end Csplus;

    record Csminus "Csminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Csminus,
      z=-1);
     annotation (preferredView = "info");
    end Csminus;

    record CsBO2 "CsBO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsBO2);
     annotation (preferredView = "info");
    end CsBO2;

    record CsBr "CsBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsBr);
     annotation (preferredView = "info");
    end CsBr;

    record CsCL "CsCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsCL);
     annotation (preferredView = "info");
    end CsCL;

    record CsF "CsF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsF);
     annotation (preferredView = "info");
    end CsF;

    record CsH "CsH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsH);
     annotation (preferredView = "info");
    end CsH;

    record CsI "CsI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsI);
     annotation (preferredView = "info");
    end CsI;

    record CsLi "CsLi(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsLi);
     annotation (preferredView = "info");
    end CsLi;

    record CsNO2 "CsNO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsNO2);
     annotation (preferredView = "info");
    end CsNO2;

    record CsNO3 "CsNO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsNO3);
     annotation (preferredView = "info");
    end CsNO3;

    record CsNa "CsNa(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsNa);
     annotation (preferredView = "info");
    end CsNa;

    record CsO "CsO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsO);
     annotation (preferredView = "info");
    end CsO;

    record CsOH "CsOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsOH);
     annotation (preferredView = "info");
    end CsOH;

    record CsRb "CsRb(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CsRb);
     annotation (preferredView = "info");
    end CsRb;

    record Cs2 "Cs2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2);
     annotation (preferredView = "info");
    end Cs2;

    record Cs2Br2 "Cs2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2Br2);
     annotation (preferredView = "info");
    end Cs2Br2;

    record Cs2CO3 "Cs2CO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2CO3);
     annotation (preferredView = "info");
    end Cs2CO3;

    record Cs2CL2 "Cs2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2CL2);
     annotation (preferredView = "info");
    end Cs2CL2;

    record Cs2F2 "Cs2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2F2);
     annotation (preferredView = "info");
    end Cs2F2;

    record Cs2I2 "Cs2I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2I2);
     annotation (preferredView = "info");
    end Cs2I2;

    record Cs2O "Cs2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2O);
     annotation (preferredView = "info");
    end Cs2O;

    record Cs2Oplus "Cs2Oplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2Oplus,
      z=1);
     annotation (preferredView = "info");
    end Cs2Oplus;

    record Cs2O2 "Cs2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2O2);
     annotation (preferredView = "info");
    end Cs2O2;

    record Cs2O2H2 "Cs2O2H2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2O2H2);
     annotation (preferredView = "info");
    end Cs2O2H2;

    record Cs2SO4 "Cs2SO4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cs2SO4);
     annotation (preferredView = "info");
    end Cs2SO4;

    record Cu "Cu(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cu);
     annotation (preferredView = "info");
    end Cu;

    record Cuplus "Cuplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cuplus,
      z=1);
     annotation (preferredView = "info");
    end Cuplus;

    record Cuminus "Cuminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cuminus,
      z=-1);
     annotation (preferredView = "info");
    end Cuminus;

    record CuCL "CuCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CuCL);
     annotation (preferredView = "info");
    end CuCL;

    record CuF "CuF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CuF);
     annotation (preferredView = "info");
    end CuF;

    record CuF2 "CuF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CuF2);
     annotation (preferredView = "info");
    end CuF2;

    record CuO "CuO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.CuO);
     annotation (preferredView = "info");
    end CuO;

    record Cu2 "Cu2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cu2);
     annotation (preferredView = "info");
    end Cu2;

    record Cu3CL3 "Cu3CL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Cu3CL3);
     annotation (preferredView = "info");
    end Cu3CL3;

    record D "D(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.D);
     annotation (preferredView = "info");
    end D;

    record Dplus "Dplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Dplus,
      z=1);
     annotation (preferredView = "info");
    end Dplus;

    record Dminus "Dminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Dminus,
      z=-1);
     annotation (preferredView = "info");
    end Dminus;

    record DBr "DBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.DBr);
     annotation (preferredView = "info");
    end DBr;

    record DCL "DCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.DCL);
     annotation (preferredView = "info");
    end DCL;

    record DF "DF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.DF);
     annotation (preferredView = "info");
    end DF;

    record DOCL "DOCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.DOCL);
     annotation (preferredView = "info");
    end DOCL;

    record DO2 "DO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.DO2);
     annotation (preferredView = "info");
    end DO2;

    record DO2minus "DO2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.DO2minus,
      z=-1);
     annotation (preferredView = "info");
    end DO2minus;

    record D2 "D2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.D2);
     annotation (preferredView = "info");
    end D2;

    record D2plus "D2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.D2plus,
      z=1);
     annotation (preferredView = "info");
    end D2plus;

    record D2minus "D2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.D2minus,
      z=-1);
     annotation (preferredView = "info");
    end D2minus;

    record D2O "D2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.D2O);
     annotation (preferredView = "info");
    end D2O;

    record D2O2 "D2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.D2O2);
     annotation (preferredView = "info");
    end D2O2;

    record D2S "D2S(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.D2S);
     annotation (preferredView = "info");
    end D2S;

    record eminus "eminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.eminus,
      z=-1);
     annotation (preferredView = "info");
    end eminus;

    record F "F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.F);
     annotation (preferredView = "info");
    end F;

    record Fplus "Fplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Fplus,
      z=1);
     annotation (preferredView = "info");
    end Fplus;

    record Fminus "Fminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Fminus,
      z=-1);
     annotation (preferredView = "info");
    end Fminus;

    record FCN "FCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FCN);
     annotation (preferredView = "info");
    end FCN;

    record FCO "FCO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FCO);
     annotation (preferredView = "info");
    end FCO;

    record FO "FO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FO);
     annotation (preferredView = "info");
    end FO;

    record FO2_FOO "FO2_FOO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FO2_FOO);
     annotation (preferredView = "info");
    end FO2_FOO;

    record FO2_OFO "FO2_OFO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FO2_OFO);
     annotation (preferredView = "info");
    end FO2_OFO;

    record F2 "F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.F2);
     annotation (preferredView = "info");
    end F2;

    record F2O "F2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.F2O);
     annotation (preferredView = "info");
    end F2O;

    record F2O2 "F2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.F2O2);
     annotation (preferredView = "info");
    end F2O2;

    record FS2F "FS2F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FS2F);
     annotation (preferredView = "info");
    end FS2F;

    record Fe "Fe(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe);
     annotation (preferredView = "info");
    end Fe;

    record Feplus "Feplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Feplus,
      z=1);
     annotation (preferredView = "info");
    end Feplus;

    record Fe_CO_5 "Fe_CO_5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe_CO_5);
     annotation (preferredView = "info");
    end Fe_CO_5;

    record FeCL "FeCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FeCL);
     annotation (preferredView = "info");
    end FeCL;

    record FeCL2 "FeCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FeCL2);
     annotation (preferredView = "info");
    end FeCL2;

    record FeCL3 "FeCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FeCL3);
     annotation (preferredView = "info");
    end FeCL3;

    record FeO "FeO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.FeO);
     annotation (preferredView = "info");
    end FeO;

    record Fe_OH_2 "Fe_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe_OH_2);
     annotation (preferredView = "info");
    end Fe_OH_2;

    record Fe2CL4 "Fe2CL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe2CL4);
     annotation (preferredView = "info");
    end Fe2CL4;

    record Fe2CL6 "Fe2CL6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Fe2CL6);
     annotation (preferredView = "info");
    end Fe2CL6;

    record Ga "Ga(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga);
     annotation (preferredView = "info");
    end Ga;

    record Gaplus "Gaplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Gaplus,
      z=1);
     annotation (preferredView = "info");
    end Gaplus;

    record GaBr "GaBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaBr);
     annotation (preferredView = "info");
    end GaBr;

    record GaBr2 "GaBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaBr2);
     annotation (preferredView = "info");
    end GaBr2;

    record GaBr3 "GaBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaBr3);
     annotation (preferredView = "info");
    end GaBr3;

    record GaCL "GaCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaCL);
     annotation (preferredView = "info");
    end GaCL;

    record GaCL2 "GaCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaCL2);
     annotation (preferredView = "info");
    end GaCL2;

    record GaCL3 "GaCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaCL3);
     annotation (preferredView = "info");
    end GaCL3;

    record GaF "GaF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaF);
     annotation (preferredView = "info");
    end GaF;

    record GaF2 "GaF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaF2);
     annotation (preferredView = "info");
    end GaF2;

    record GaF3 "GaF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaF3);
     annotation (preferredView = "info");
    end GaF3;

    record GaH "GaH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaH);
     annotation (preferredView = "info");
    end GaH;

    record GaI "GaI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaI);
     annotation (preferredView = "info");
    end GaI;

    record GaI2 "GaI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaI2);
     annotation (preferredView = "info");
    end GaI2;

    record GaI3 "GaI3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaI3);
     annotation (preferredView = "info");
    end GaI3;

    record GaO "GaO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaO);
     annotation (preferredView = "info");
    end GaO;

    record GaOH "GaOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GaOH);
     annotation (preferredView = "info");
    end GaOH;

    record Ga2Br2 "Ga2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2Br2);
     annotation (preferredView = "info");
    end Ga2Br2;

    record Ga2Br4 "Ga2Br4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2Br4);
     annotation (preferredView = "info");
    end Ga2Br4;

    record Ga2Br6 "Ga2Br6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2Br6);
     annotation (preferredView = "info");
    end Ga2Br6;

    record Ga2CL2 "Ga2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2CL2);
     annotation (preferredView = "info");
    end Ga2CL2;

    record Ga2CL4 "Ga2CL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2CL4);
     annotation (preferredView = "info");
    end Ga2CL4;

    record Ga2CL6 "Ga2CL6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2CL6);
     annotation (preferredView = "info");
    end Ga2CL6;

    record Ga2F2 "Ga2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2F2);
     annotation (preferredView = "info");
    end Ga2F2;

    record Ga2F4 "Ga2F4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2F4);
     annotation (preferredView = "info");
    end Ga2F4;

    record Ga2F6 "Ga2F6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2F6);
     annotation (preferredView = "info");
    end Ga2F6;

    record Ga2I2 "Ga2I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2I2);
     annotation (preferredView = "info");
    end Ga2I2;

    record Ga2I4 "Ga2I4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2I4);
     annotation (preferredView = "info");
    end Ga2I4;

    record Ga2I6 "Ga2I6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2I6);
     annotation (preferredView = "info");
    end Ga2I6;

    record Ga2O "Ga2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ga2O);
     annotation (preferredView = "info");
    end Ga2O;

    record Ge "Ge(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ge);
     annotation (preferredView = "info");
    end Ge;

    record Geplus "Geplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Geplus,
      z=1);
     annotation (preferredView = "info");
    end Geplus;

    record Geminus "Geminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Geminus,
      z=-1);
     annotation (preferredView = "info");
    end Geminus;

    record GeBr "GeBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeBr);
     annotation (preferredView = "info");
    end GeBr;

    record GeBr2 "GeBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeBr2);
     annotation (preferredView = "info");
    end GeBr2;

    record GeBr3 "GeBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeBr3);
     annotation (preferredView = "info");
    end GeBr3;

    record GeBr4 "GeBr4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeBr4);
     annotation (preferredView = "info");
    end GeBr4;

    record GeCL "GeCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeCL);
     annotation (preferredView = "info");
    end GeCL;

    record GeCL2 "GeCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeCL2);
     annotation (preferredView = "info");
    end GeCL2;

    record GeCL3 "GeCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeCL3);
     annotation (preferredView = "info");
    end GeCL3;

    record GeCL4 "GeCL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeCL4);
     annotation (preferredView = "info");
    end GeCL4;

    record GeF "GeF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeF);
     annotation (preferredView = "info");
    end GeF;

    record GeF2 "GeF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeF2);
     annotation (preferredView = "info");
    end GeF2;

    record GeF3 "GeF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeF3);
     annotation (preferredView = "info");
    end GeF3;

    record GeF4 "GeF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeF4);
     annotation (preferredView = "info");
    end GeF4;

    record GeH4 "GeH4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeH4);
     annotation (preferredView = "info");
    end GeH4;

    record GeI "GeI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeI);
     annotation (preferredView = "info");
    end GeI;

    record GeO "GeO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeO);
     annotation (preferredView = "info");
    end GeO;

    record GeO2 "GeO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeO2);
     annotation (preferredView = "info");
    end GeO2;

    record GeS "GeS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeS);
     annotation (preferredView = "info");
    end GeS;

    record GeS2 "GeS2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.GeS2);
     annotation (preferredView = "info");
    end GeS2;

    record Ge2 "Ge2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ge2);
     annotation (preferredView = "info");
    end Ge2;

    record H "H(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H);
     annotation (preferredView = "info");
    end H;

    record Hplus "Hplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Hplus,
      z=1);
     annotation (preferredView = "info");
    end Hplus;

    record Hminus "Hminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Hminus,
      z=-1);
     annotation (preferredView = "info");
    end Hminus;

    record HALO "HALO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HALO);
     annotation (preferredView = "info");
    end HALO;

    record HALO2 "HALO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HALO2);
     annotation (preferredView = "info");
    end HALO2;

    record HBO "HBO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HBO);
     annotation (preferredView = "info");
    end HBO;

    record HBOplus "HBOplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HBOplus,
      z=1);
     annotation (preferredView = "info");
    end HBOplus;

    record HBO2 "HBO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HBO2);
     annotation (preferredView = "info");
    end HBO2;

    record HBS "HBS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HBS);
     annotation (preferredView = "info");
    end HBS;

    record HBSplus "HBSplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HBSplus,
      z=1);
     annotation (preferredView = "info");
    end HBSplus;

    record HCN "HCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HCN);
     annotation (preferredView = "info");
    end HCN;

    record HCO "HCO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HCO);
     annotation (preferredView = "info");
    end HCO;

    record HCOplus "HCOplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HCOplus,
      z=1);
     annotation (preferredView = "info");
    end HCOplus;

    record HCCN "HCCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HCCN);
     annotation (preferredView = "info");
    end HCCN;

    record HCCO "HCCO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HCCO);
     annotation (preferredView = "info");
    end HCCO;

    record HCL "HCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HCL);
     annotation (preferredView = "info");
    end HCL;

    record HD "HD(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HD);
     annotation (preferredView = "info");
    end HD;

    record HDplus "HDplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HDplus,
      z=1);
     annotation (preferredView = "info");
    end HDplus;

    record HDO "HDO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HDO);
     annotation (preferredView = "info");
    end HDO;

    record HDO2 "HDO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HDO2);
     annotation (preferredView = "info");
    end HDO2;

    record HF "HF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HF);
     annotation (preferredView = "info");
    end HF;

    record HI "HI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HI);
     annotation (preferredView = "info");
    end HI;

    record HNC "HNC(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HNC);
     annotation (preferredView = "info");
    end HNC;

    record HNCO "HNCO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HNCO);
     annotation (preferredView = "info");
    end HNCO;

    record HNO "HNO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HNO);
     annotation (preferredView = "info");
    end HNO;

    record HNO2 "HNO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HNO2);
     annotation (preferredView = "info");
    end HNO2;

    record HNO3 "HNO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HNO3);
     annotation (preferredView = "info");
    end HNO3;

    record HOCL "HOCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HOCL);
     annotation (preferredView = "info");
    end HOCL;

    record HOF "HOF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HOF);
     annotation (preferredView = "info");
    end HOF;

    record HO2 "HO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HO2);
     annotation (preferredView = "info");
    end HO2;

    record HO2minus "HO2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HO2minus,
      z=-1);
     annotation (preferredView = "info");
    end HO2minus;

    record HPO "HPO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HPO);
     annotation (preferredView = "info");
    end HPO;

    record HSO3F "HSO3F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HSO3F);
     annotation (preferredView = "info");
    end HSO3F;

    record H2 "H2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2);
     annotation (preferredView = "info");
    end H2;

    record H2plus "H2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2plus,
      z=1);
     annotation (preferredView = "info");
    end H2plus;

    record H2minus "H2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2minus,
      z=-1);
     annotation (preferredView = "info");
    end H2minus;

    record HBOH "HBOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HBOH);
     annotation (preferredView = "info");
    end HBOH;

    record HCOOH "HCOOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HCOOH);
     annotation (preferredView = "info");
    end HCOOH;

    record H2F2 "H2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2F2);
     annotation (preferredView = "info");
    end H2F2;

    record H2O "H2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2O);
     annotation (preferredView = "info");
    end H2O;

    record H2Oplus "H2Oplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2Oplus,
      z=1);
     annotation (preferredView = "info");
    end H2Oplus;

    record H2O2 "H2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2O2);
     annotation (preferredView = "info");
    end H2O2;

    record H2S "H2S(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2S);
     annotation (preferredView = "info");
    end H2S;

    record H2SO4 "H2SO4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2SO4);
     annotation (preferredView = "info");
    end H2SO4;

    record H2BOH "H2BOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H2BOH);
     annotation (preferredView = "info");
    end H2BOH;

    record HB_OH_2 "HB_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HB_OH_2);
     annotation (preferredView = "info");
    end HB_OH_2;

    record H3BO3 "H3BO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H3BO3);
     annotation (preferredView = "info");
    end H3BO3;

    record H3B3O3 "H3B3O3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H3B3O3);
     annotation (preferredView = "info");
    end H3B3O3;

    record H3B3O6 "H3B3O6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H3B3O6);
     annotation (preferredView = "info");
    end H3B3O6;

    record H3F3 "H3F3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H3F3);
     annotation (preferredView = "info");
    end H3F3;

    record H3Oplus "H3Oplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H3Oplus,
      z=1);
     annotation (preferredView = "info");
    end H3Oplus;

    record H4F4 "H4F4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H4F4);
     annotation (preferredView = "info");
    end H4F4;

    record H5F5 "H5F5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H5F5);
     annotation (preferredView = "info");
    end H5F5;

    record H6F6 "H6F6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H6F6);
     annotation (preferredView = "info");
    end H6F6;

    record H7F7 "H7F7(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.H7F7);
     annotation (preferredView = "info");
    end H7F7;

    record He "He(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.He);
     annotation (preferredView = "info");
    end He;

    record Heplus "Heplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Heplus,
      z=1);
     annotation (preferredView = "info");
    end Heplus;

    record Hg "Hg(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Hg);
     annotation (preferredView = "info");
    end Hg;

    record Hgplus "Hgplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Hgplus,
      z=1);
     annotation (preferredView = "info");
    end Hgplus;

    record HgBr2 "HgBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.HgBr2);
     annotation (preferredView = "info");
    end HgBr2;

    record I "I(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.I);
     annotation (preferredView = "info");
    end I;

    record Iplus "Iplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Iplus,
      z=1);
     annotation (preferredView = "info");
    end Iplus;

    record Iminus "Iminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Iminus,
      z=-1);
     annotation (preferredView = "info");
    end Iminus;

    record IF5 "IF5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.IF5);
     annotation (preferredView = "info");
    end IF5;

    record IF7 "IF7(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.IF7);
     annotation (preferredView = "info");
    end IF7;

    record I2 "I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.I2);
     annotation (preferredView = "info");
    end I2;

    record In "In(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In);
     annotation (preferredView = "info");
    end In;

    record Inplus "Inplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Inplus,
      z=1);
     annotation (preferredView = "info");
    end Inplus;

    record InBr "InBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InBr);
     annotation (preferredView = "info");
    end InBr;

    record InBr2 "InBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InBr2);
     annotation (preferredView = "info");
    end InBr2;

    record InBr3 "InBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InBr3);
     annotation (preferredView = "info");
    end InBr3;

    record InCL "InCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InCL);
     annotation (preferredView = "info");
    end InCL;

    record InCL2 "InCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InCL2);
     annotation (preferredView = "info");
    end InCL2;

    record InCL3 "InCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InCL3);
     annotation (preferredView = "info");
    end InCL3;

    record InF "InF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InF);
     annotation (preferredView = "info");
    end InF;

    record InF2 "InF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InF2);
     annotation (preferredView = "info");
    end InF2;

    record InF3 "InF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InF3);
     annotation (preferredView = "info");
    end InF3;

    record InH "InH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InH);
     annotation (preferredView = "info");
    end InH;

    record InI "InI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InI);
     annotation (preferredView = "info");
    end InI;

    record InI2 "InI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InI2);
     annotation (preferredView = "info");
    end InI2;

    record InI3 "InI3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InI3);
     annotation (preferredView = "info");
    end InI3;

    record InO "InO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InO);
     annotation (preferredView = "info");
    end InO;

    record InOH "InOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.InOH);
     annotation (preferredView = "info");
    end InOH;

    record In2Br2 "In2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2Br2);
     annotation (preferredView = "info");
    end In2Br2;

    record In2Br4 "In2Br4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2Br4);
     annotation (preferredView = "info");
    end In2Br4;

    record In2Br6 "In2Br6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2Br6);
     annotation (preferredView = "info");
    end In2Br6;

    record In2CL2 "In2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2CL2);
     annotation (preferredView = "info");
    end In2CL2;

    record In2CL4 "In2CL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2CL4);
     annotation (preferredView = "info");
    end In2CL4;

    record In2CL6 "In2CL6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2CL6);
     annotation (preferredView = "info");
    end In2CL6;

    record In2F2 "In2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2F2);
     annotation (preferredView = "info");
    end In2F2;

    record In2F4 "In2F4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2F4);
     annotation (preferredView = "info");
    end In2F4;

    record In2F6 "In2F6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2F6);
     annotation (preferredView = "info");
    end In2F6;

    record In2I2 "In2I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2I2);
     annotation (preferredView = "info");
    end In2I2;

    record In2I4 "In2I4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2I4);
     annotation (preferredView = "info");
    end In2I4;

    record In2I6 "In2I6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2I6);
     annotation (preferredView = "info");
    end In2I6;

    record In2O "In2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.In2O);
     annotation (preferredView = "info");
    end In2O;

    record K "K(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K);
     annotation (preferredView = "info");
    end K;

    record Kplus "Kplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Kplus,
      z=1);
     annotation (preferredView = "info");
    end Kplus;

    record Kminus "Kminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Kminus,
      z=-1);
     annotation (preferredView = "info");
    end Kminus;

    record KALF4 "KALF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KALF4);
     annotation (preferredView = "info");
    end KALF4;

    record KBO2 "KBO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KBO2);
     annotation (preferredView = "info");
    end KBO2;

    record KBr "KBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KBr);
     annotation (preferredView = "info");
    end KBr;

    record KCN "KCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KCN);
     annotation (preferredView = "info");
    end KCN;

    record KCL "KCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KCL);
     annotation (preferredView = "info");
    end KCL;

    record KF "KF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KF);
     annotation (preferredView = "info");
    end KF;

    record KH "KH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KH);
     annotation (preferredView = "info");
    end KH;

    record KI "KI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KI);
     annotation (preferredView = "info");
    end KI;

    record KLi "KLi(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KLi);
     annotation (preferredView = "info");
    end KLi;

    record KNO2 "KNO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KNO2);
     annotation (preferredView = "info");
    end KNO2;

    record KNO3 "KNO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KNO3);
     annotation (preferredView = "info");
    end KNO3;

    record KNa "KNa(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KNa);
     annotation (preferredView = "info");
    end KNa;

    record KO "KO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KO);
     annotation (preferredView = "info");
    end KO;

    record KOH "KOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.KOH);
     annotation (preferredView = "info");
    end KOH;

    record K2 "K2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2);
     annotation (preferredView = "info");
    end K2;

    record K2plus "K2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2plus,
      z=1);
     annotation (preferredView = "info");
    end K2plus;

    record K2Br2 "K2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2Br2);
     annotation (preferredView = "info");
    end K2Br2;

    record K2CO3 "K2CO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2CO3);
     annotation (preferredView = "info");
    end K2CO3;

    record K2C2N2 "K2C2N2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2C2N2);
     annotation (preferredView = "info");
    end K2C2N2;

    record K2CL2 "K2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2CL2);
     annotation (preferredView = "info");
    end K2CL2;

    record K2F2 "K2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2F2);
     annotation (preferredView = "info");
    end K2F2;

    record K2I2 "K2I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2I2);
     annotation (preferredView = "info");
    end K2I2;

    record K2O "K2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2O);
     annotation (preferredView = "info");
    end K2O;

    record K2Oplus "K2Oplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2Oplus,
      z=1);
     annotation (preferredView = "info");
    end K2Oplus;

    record K2O2 "K2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2O2);
     annotation (preferredView = "info");
    end K2O2;

    record K2O2H2 "K2O2H2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2O2H2);
     annotation (preferredView = "info");
    end K2O2H2;

    record K2SO4 "K2SO4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.K2SO4);
     annotation (preferredView = "info");
    end K2SO4;

    record Kr "Kr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Kr);
     annotation (preferredView = "info");
    end Kr;

    record Krplus "Krplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Krplus,
      z=1);
     annotation (preferredView = "info");
    end Krplus;

    record Li "Li(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li);
     annotation (preferredView = "info");
    end Li;

    record Liplus "Liplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Liplus,
      z=1);
     annotation (preferredView = "info");
    end Liplus;

    record Liminus "Liminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Liminus,
      z=-1);
     annotation (preferredView = "info");
    end Liminus;

    record LiALF4 "LiALF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiALF4);
     annotation (preferredView = "info");
    end LiALF4;

    record LiBO2 "LiBO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiBO2);
     annotation (preferredView = "info");
    end LiBO2;

    record LiBr "LiBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiBr);
     annotation (preferredView = "info");
    end LiBr;

    record LiCL "LiCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiCL);
     annotation (preferredView = "info");
    end LiCL;

    record LiF "LiF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiF);
     annotation (preferredView = "info");
    end LiF;

    record LiH "LiH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiH);
     annotation (preferredView = "info");
    end LiH;

    record LiI "LiI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiI);
     annotation (preferredView = "info");
    end LiI;

    record LiN "LiN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiN);
     annotation (preferredView = "info");
    end LiN;

    record LiNO2 "LiNO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiNO2);
     annotation (preferredView = "info");
    end LiNO2;

    record LiNO3 "LiNO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiNO3);
     annotation (preferredView = "info");
    end LiNO3;

    record LiO "LiO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiO);
     annotation (preferredView = "info");
    end LiO;

    record LiOF "LiOF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiOF);
     annotation (preferredView = "info");
    end LiOF;

    record LiOH "LiOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiOH);
     annotation (preferredView = "info");
    end LiOH;

    record LiON "LiON(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.LiON);
     annotation (preferredView = "info");
    end LiON;

    record Li2 "Li2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2);
     annotation (preferredView = "info");
    end Li2;

    record Li2plus "Li2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2plus,
      z=1);
     annotation (preferredView = "info");
    end Li2plus;

    record Li2Br2 "Li2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2Br2);
     annotation (preferredView = "info");
    end Li2Br2;

    record Li2F2 "Li2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2F2);
     annotation (preferredView = "info");
    end Li2F2;

    record Li2I2 "Li2I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2I2);
     annotation (preferredView = "info");
    end Li2I2;

    record Li2O "Li2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2O);
     annotation (preferredView = "info");
    end Li2O;

    record Li2Oplus "Li2Oplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2Oplus,
      z=1);
     annotation (preferredView = "info");
    end Li2Oplus;

    record Li2O2 "Li2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2O2);
     annotation (preferredView = "info");
    end Li2O2;

    record Li2O2H2 "Li2O2H2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2O2H2);
     annotation (preferredView = "info");
    end Li2O2H2;

    record Li2SO4 "Li2SO4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li2SO4);
     annotation (preferredView = "info");
    end Li2SO4;

    record Li3plus "Li3plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3plus,
      z=1);
     annotation (preferredView = "info");
    end Li3plus;

    record Li3Br3 "Li3Br3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3Br3);
     annotation (preferredView = "info");
    end Li3Br3;

    record Li3CL3 "Li3CL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3CL3);
     annotation (preferredView = "info");
    end Li3CL3;

    record Li3F3 "Li3F3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3F3);
     annotation (preferredView = "info");
    end Li3F3;

    record Li3I3 "Li3I3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Li3I3);
     annotation (preferredView = "info");
    end Li3I3;

    record Mg "Mg(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mg);
     annotation (preferredView = "info");
    end Mg;

    record Mgplus "Mgplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mgplus,
      z=1);
     annotation (preferredView = "info");
    end Mgplus;

    record MgBr "MgBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgBr);
     annotation (preferredView = "info");
    end MgBr;

    record MgBr2 "MgBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgBr2);
     annotation (preferredView = "info");
    end MgBr2;

    record MgCL "MgCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgCL);
     annotation (preferredView = "info");
    end MgCL;

    record MgCLplus "MgCLplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgCLplus,
      z=1);
     annotation (preferredView = "info");
    end MgCLplus;

    record MgCL2 "MgCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgCL2);
     annotation (preferredView = "info");
    end MgCL2;

    record MgF "MgF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgF);
     annotation (preferredView = "info");
    end MgF;

    record MgFplus "MgFplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgFplus,
      z=1);
     annotation (preferredView = "info");
    end MgFplus;

    record MgF2 "MgF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgF2);
     annotation (preferredView = "info");
    end MgF2;

    record MgF2plus "MgF2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgF2plus,
      z=1);
     annotation (preferredView = "info");
    end MgF2plus;

    record MgH "MgH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgH);
     annotation (preferredView = "info");
    end MgH;

    record MgI "MgI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgI);
     annotation (preferredView = "info");
    end MgI;

    record MgI2 "MgI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgI2);
     annotation (preferredView = "info");
    end MgI2;

    record MgN "MgN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgN);
     annotation (preferredView = "info");
    end MgN;

    record MgO "MgO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgO);
     annotation (preferredView = "info");
    end MgO;

    record MgOH "MgOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgOH);
     annotation (preferredView = "info");
    end MgOH;

    record MgOHplus "MgOHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgOHplus,
      z=1);
     annotation (preferredView = "info");
    end MgOHplus;

    record Mg_OH_2 "Mg_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mg_OH_2);
     annotation (preferredView = "info");
    end Mg_OH_2;

    record MgS "MgS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MgS);
     annotation (preferredView = "info");
    end MgS;

    record Mg2 "Mg2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mg2);
     annotation (preferredView = "info");
    end Mg2;

    record Mg2F4 "Mg2F4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mg2F4);
     annotation (preferredView = "info");
    end Mg2F4;

    record Mn "Mn(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mn);
     annotation (preferredView = "info");
    end Mn;

    record Mnplus "Mnplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mnplus,
      z=1);
     annotation (preferredView = "info");
    end Mnplus;

    record Mo "Mo(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo);
     annotation (preferredView = "info");
    end Mo;

    record Moplus "Moplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Moplus,
      z=1);
     annotation (preferredView = "info");
    end Moplus;

    record Mominus "Mominus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mominus,
      z=-1);
     annotation (preferredView = "info");
    end Mominus;

    record MoO "MoO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MoO);
     annotation (preferredView = "info");
    end MoO;

    record MoO2 "MoO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MoO2);
     annotation (preferredView = "info");
    end MoO2;

    record MoO3 "MoO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MoO3);
     annotation (preferredView = "info");
    end MoO3;

    record MoO3minus "MoO3minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.MoO3minus,
      z=-1);
     annotation (preferredView = "info");
    end MoO3minus;

    record Mo2O6 "Mo2O6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo2O6);
     annotation (preferredView = "info");
    end Mo2O6;

    record Mo3O9 "Mo3O9(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo3O9);
     annotation (preferredView = "info");
    end Mo3O9;

    record Mo4O12 "Mo4O12(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo4O12);
     annotation (preferredView = "info");
    end Mo4O12;

    record Mo5O15 "Mo5O15(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Mo5O15);
     annotation (preferredView = "info");
    end Mo5O15;

    record N "N(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N);
     annotation (preferredView = "info");
    end N;

    record Nplus "Nplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Nplus,
      z=1);
     annotation (preferredView = "info");
    end Nplus;

    record Nminus "Nminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Nminus,
      z=-1);
     annotation (preferredView = "info");
    end Nminus;

    record NCO "NCO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NCO);
     annotation (preferredView = "info");
    end NCO;

    record ND "ND(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ND);
     annotation (preferredView = "info");
    end ND;

    record ND2 "ND2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ND2);
     annotation (preferredView = "info");
    end ND2;

    record ND3 "ND3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ND3);
     annotation (preferredView = "info");
    end ND3;

    record NF "NF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NF);
     annotation (preferredView = "info");
    end NF;

    record NF2 "NF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NF2);
     annotation (preferredView = "info");
    end NF2;

    record NF3 "NF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NF3);
     annotation (preferredView = "info");
    end NF3;

    record NH "NH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NH);
     annotation (preferredView = "info");
    end NH;

    record NHplus "NHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NHplus,
      z=1);
     annotation (preferredView = "info");
    end NHplus;

    record NHF "NHF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NHF);
     annotation (preferredView = "info");
    end NHF;

    record NHF2 "NHF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NHF2);
     annotation (preferredView = "info");
    end NHF2;

    record NH2 "NH2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NH2);
     annotation (preferredView = "info");
    end NH2;

    record NH2F "NH2F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NH2F);
     annotation (preferredView = "info");
    end NH2F;

    record NH3 "NH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NH3);
     annotation (preferredView = "info");
    end NH3;

    record NH2OH "NH2OH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NH2OH);
     annotation (preferredView = "info");
    end NH2OH;

    record NH4plus "NH4plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NH4plus,
      z=1);
     annotation (preferredView = "info");
    end NH4plus;

    record NO "NO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NO);
     annotation (preferredView = "info");
    end NO;

    record NOCL "NOCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NOCL);
     annotation (preferredView = "info");
    end NOCL;

    record NOF "NOF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NOF);
     annotation (preferredView = "info");
    end NOF;

    record NOF3 "NOF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NOF3);
     annotation (preferredView = "info");
    end NOF3;

    record NO2 "NO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NO2);
     annotation (preferredView = "info");
    end NO2;

    record NO2minus "NO2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NO2minus,
      z=-1);
     annotation (preferredView = "info");
    end NO2minus;

    record NO2CL "NO2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NO2CL);
     annotation (preferredView = "info");
    end NO2CL;

    record NO2F "NO2F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NO2F);
     annotation (preferredView = "info");
    end NO2F;

    record NO3 "NO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NO3);
     annotation (preferredView = "info");
    end NO3;

    record NO3minus "NO3minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NO3minus,
      z=-1);
     annotation (preferredView = "info");
    end NO3minus;

    record NO3F "NO3F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NO3F);
     annotation (preferredView = "info");
    end NO3F;

    record N2 "N2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2);
     annotation (preferredView = "info");
    end N2;

    record N2plus "N2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2plus,
      z=1);
     annotation (preferredView = "info");
    end N2plus;

    record N2minus "N2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2minus,
      z=-1);
     annotation (preferredView = "info");
    end N2minus;

    record NCN "NCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NCN);
     annotation (preferredView = "info");
    end NCN;

    record N2D2_cis "N2D2_cis(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2D2_cis);
     annotation (preferredView = "info");
    end N2D2_cis;

    record N2F2 "N2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2F2);
     annotation (preferredView = "info");
    end N2F2;

    record N2F4 "N2F4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2F4);
     annotation (preferredView = "info");
    end N2F4;

    record N2H2 "N2H2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2H2);
     annotation (preferredView = "info");
    end N2H2;

    record NH2NO2 "NH2NO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NH2NO2);
     annotation (preferredView = "info");
    end NH2NO2;

    record N2H4 "N2H4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2H4);
     annotation (preferredView = "info");
    end N2H4;

    record N2O "N2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2O);
     annotation (preferredView = "info");
    end N2O;

    record N2Oplus "N2Oplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2Oplus,
      z=1);
     annotation (preferredView = "info");
    end N2Oplus;

    record N2O3 "N2O3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2O3);
     annotation (preferredView = "info");
    end N2O3;

    record N2O4 "N2O4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2O4);
     annotation (preferredView = "info");
    end N2O4;

    record N2O5 "N2O5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N2O5);
     annotation (preferredView = "info");
    end N2O5;

    record N3 "N3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N3);
     annotation (preferredView = "info");
    end N3;

    record N3H "N3H(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.N3H);
     annotation (preferredView = "info");
    end N3H;

    record Na "Na(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na);
     annotation (preferredView = "info");
    end Na;

    record Naplus "Naplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Naplus,
      z=1);
     annotation (preferredView = "info");
    end Naplus;

    record Naminus "Naminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Naminus,
      z=-1);
     annotation (preferredView = "info");
    end Naminus;

    record NaALF4 "NaALF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaALF4);
     annotation (preferredView = "info");
    end NaALF4;

    record NaBO2 "NaBO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaBO2);
     annotation (preferredView = "info");
    end NaBO2;

    record NaBr "NaBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaBr);
     annotation (preferredView = "info");
    end NaBr;

    record NaCN "NaCN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaCN);
     annotation (preferredView = "info");
    end NaCN;

    record NaCL "NaCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaCL);
     annotation (preferredView = "info");
    end NaCL;

    record NaF "NaF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaF);
     annotation (preferredView = "info");
    end NaF;

    record NaH "NaH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaH);
     annotation (preferredView = "info");
    end NaH;

    record NaI "NaI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaI);
     annotation (preferredView = "info");
    end NaI;

    record NaLi "NaLi(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaLi);
     annotation (preferredView = "info");
    end NaLi;

    record NaNO2 "NaNO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaNO2);
     annotation (preferredView = "info");
    end NaNO2;

    record NaNO3 "NaNO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaNO3);
     annotation (preferredView = "info");
    end NaNO3;

    record NaO "NaO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaO);
     annotation (preferredView = "info");
    end NaO;

    record NaOH "NaOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaOH);
     annotation (preferredView = "info");
    end NaOH;

    record NaOHplus "NaOHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NaOHplus,
      z=1);
     annotation (preferredView = "info");
    end NaOHplus;

    record Na2 "Na2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2);
     annotation (preferredView = "info");
    end Na2;

    record Na2Br2 "Na2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2Br2);
     annotation (preferredView = "info");
    end Na2Br2;

    record Na2CL2 "Na2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2CL2);
     annotation (preferredView = "info");
    end Na2CL2;

    record Na2F2 "Na2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2F2);
     annotation (preferredView = "info");
    end Na2F2;

    record Na2I2 "Na2I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2I2);
     annotation (preferredView = "info");
    end Na2I2;

    record Na2O "Na2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2O);
     annotation (preferredView = "info");
    end Na2O;

    record Na2Oplus "Na2Oplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2Oplus,
      z=1);
     annotation (preferredView = "info");
    end Na2Oplus;

    record Na2O2 "Na2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2O2);
     annotation (preferredView = "info");
    end Na2O2;

    record Na2O2H2 "Na2O2H2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2O2H2);
     annotation (preferredView = "info");
    end Na2O2H2;

    record Na2SO4 "Na2SO4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na2SO4);
     annotation (preferredView = "info");
    end Na2SO4;

    record Na3CL3 "Na3CL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na3CL3);
     annotation (preferredView = "info");
    end Na3CL3;

    record Na3F3 "Na3F3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Na3F3);
     annotation (preferredView = "info");
    end Na3F3;

    record Nb "Nb(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Nb);
     annotation (preferredView = "info");
    end Nb;

    record Nbplus "Nbplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Nbplus,
      z=1);
     annotation (preferredView = "info");
    end Nbplus;

    record Nbminus "Nbminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Nbminus,
      z=-1);
     annotation (preferredView = "info");
    end Nbminus;

    record NbCL5 "NbCL5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NbCL5);
     annotation (preferredView = "info");
    end NbCL5;

    record NbO "NbO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NbO);
     annotation (preferredView = "info");
    end NbO;

    record NbOCL3 "NbOCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NbOCL3);
     annotation (preferredView = "info");
    end NbOCL3;

    record NbO2 "NbO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NbO2);
     annotation (preferredView = "info");
    end NbO2;

    record Ne "Ne(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ne);
     annotation (preferredView = "info");
    end Ne;

    record Neplus "Neplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Neplus,
      z=1);
     annotation (preferredView = "info");
    end Neplus;

    record Ni "Ni(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ni);
     annotation (preferredView = "info");
    end Ni;

    record Niplus "Niplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Niplus,
      z=1);
     annotation (preferredView = "info");
    end Niplus;

    record Niminus "Niminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Niminus,
      z=-1);
     annotation (preferredView = "info");
    end Niminus;

    record NiCL "NiCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NiCL);
     annotation (preferredView = "info");
    end NiCL;

    record NiCL2 "NiCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NiCL2);
     annotation (preferredView = "info");
    end NiCL2;

    record NiO "NiO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NiO);
     annotation (preferredView = "info");
    end NiO;

    record NiS "NiS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.NiS);
     annotation (preferredView = "info");
    end NiS;

    record O "O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.O);
     annotation (preferredView = "info");
    end O;

    record Oplus "Oplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Oplus,
      z=1);
     annotation (preferredView = "info");
    end Oplus;

    record Ominus "Ominus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ominus,
      z=-1);
     annotation (preferredView = "info");
    end Ominus;

    record OD "OD(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.OD);
     annotation (preferredView = "info");
    end OD;

    record ODminus "ODminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ODminus,
      z=-1);
     annotation (preferredView = "info");
    end ODminus;

    record OH "OH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.OH);
     annotation (preferredView = "info");
    end OH;

    record OHplus "OHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.OHplus,
      z=1);
     annotation (preferredView = "info");
    end OHplus;

    record OHminus "OHminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.OHminus,
      z=-1);
     annotation (preferredView = "info");
    end OHminus;

    record O2 "O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.O2);
     annotation (preferredView = "info");
    end O2;

    record O2plus "O2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.O2plus,
      z=1);
     annotation (preferredView = "info");
    end O2plus;

    record O2minus "O2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.O2minus,
      z=-1);
     annotation (preferredView = "info");
    end O2minus;

    record O3 "O3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.O3);
     annotation (preferredView = "info");
    end O3;

    record P "P(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P);
     annotation (preferredView = "info");
    end P;

    record Pplus "Pplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Pplus,
      z=1);
     annotation (preferredView = "info");
    end Pplus;

    record Pminus "Pminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Pminus,
      z=-1);
     annotation (preferredView = "info");
    end Pminus;

    record PCL "PCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL);
     annotation (preferredView = "info");
    end PCL;

    record PCL2 "PCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL2);
     annotation (preferredView = "info");
    end PCL2;

    record PCL2minus "PCL2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL2minus,
      z=-1);
     annotation (preferredView = "info");
    end PCL2minus;

    record PCL3 "PCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL3);
     annotation (preferredView = "info");
    end PCL3;

    record PCL5 "PCL5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PCL5);
     annotation (preferredView = "info");
    end PCL5;

    record PF "PF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PF);
     annotation (preferredView = "info");
    end PF;

    record PFplus "PFplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PFplus,
      z=1);
     annotation (preferredView = "info");
    end PFplus;

    record PFminus "PFminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PFminus,
      z=-1);
     annotation (preferredView = "info");
    end PFminus;

    record PFCL "PFCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PFCL);
     annotation (preferredView = "info");
    end PFCL;

    record PFCLminus "PFCLminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PFCLminus,
      z=-1);
     annotation (preferredView = "info");
    end PFCLminus;

    record PFCL2 "PFCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PFCL2);
     annotation (preferredView = "info");
    end PFCL2;

    record PFCL4 "PFCL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PFCL4);
     annotation (preferredView = "info");
    end PFCL4;

    record PF2 "PF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PF2);
     annotation (preferredView = "info");
    end PF2;

    record PF2minus "PF2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PF2minus,
      z=-1);
     annotation (preferredView = "info");
    end PF2minus;

    record PF2CL "PF2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PF2CL);
     annotation (preferredView = "info");
    end PF2CL;

    record PF2CL3 "PF2CL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PF2CL3);
     annotation (preferredView = "info");
    end PF2CL3;

    record PF3 "PF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PF3);
     annotation (preferredView = "info");
    end PF3;

    record PF3CL2 "PF3CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PF3CL2);
     annotation (preferredView = "info");
    end PF3CL2;

    record PF4CL "PF4CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PF4CL);
     annotation (preferredView = "info");
    end PF4CL;

    record PF5 "PF5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PF5);
     annotation (preferredView = "info");
    end PF5;

    record PH "PH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PH);
     annotation (preferredView = "info");
    end PH;

    record PH2 "PH2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PH2);
     annotation (preferredView = "info");
    end PH2;

    record PH2minus "PH2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PH2minus,
      z=-1);
     annotation (preferredView = "info");
    end PH2minus;

    record PH3 "PH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PH3);
     annotation (preferredView = "info");
    end PH3;

    record PN "PN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PN);
     annotation (preferredView = "info");
    end PN;

    record PO "PO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PO);
     annotation (preferredView = "info");
    end PO;

    record POminus "POminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.POminus,
      z=-1);
     annotation (preferredView = "info");
    end POminus;

    record POCL3 "POCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.POCL3);
     annotation (preferredView = "info");
    end POCL3;

    record POFCL2 "POFCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.POFCL2);
     annotation (preferredView = "info");
    end POFCL2;

    record POF2CL "POF2CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.POF2CL);
     annotation (preferredView = "info");
    end POF2CL;

    record POF3 "POF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.POF3);
     annotation (preferredView = "info");
    end POF3;

    record PO2 "PO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PO2);
     annotation (preferredView = "info");
    end PO2;

    record PO2minus "PO2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PO2minus,
      z=-1);
     annotation (preferredView = "info");
    end PO2minus;

    record PS "PS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PS);
     annotation (preferredView = "info");
    end PS;

    record P2 "P2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P2);
     annotation (preferredView = "info");
    end P2;

    record P2O3 "P2O3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P2O3);
     annotation (preferredView = "info");
    end P2O3;

    record P2O4 "P2O4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P2O4);
     annotation (preferredView = "info");
    end P2O4;

    record P2O5 "P2O5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P2O5);
     annotation (preferredView = "info");
    end P2O5;

    record P3 "P3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P3);
     annotation (preferredView = "info");
    end P3;

    record P3O6 "P3O6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P3O6);
     annotation (preferredView = "info");
    end P3O6;

    record P4 "P4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P4);
     annotation (preferredView = "info");
    end P4;

    record P4O6 "P4O6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O6);
     annotation (preferredView = "info");
    end P4O6;

    record P4O7 "P4O7(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O7);
     annotation (preferredView = "info");
    end P4O7;

    record P4O8 "P4O8(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O8);
     annotation (preferredView = "info");
    end P4O8;

    record P4O9 "P4O9(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O9);
     annotation (preferredView = "info");
    end P4O9;

    record P4O10 "P4O10(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.P4O10);
     annotation (preferredView = "info");
    end P4O10;

    record Pb "Pb(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Pb);
     annotation (preferredView = "info");
    end Pb;

    record Pbplus "Pbplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Pbplus,
      z=1);
     annotation (preferredView = "info");
    end Pbplus;

    record Pbminus "Pbminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Pbminus,
      z=-1);
     annotation (preferredView = "info");
    end Pbminus;

    record PbBr "PbBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbBr);
     annotation (preferredView = "info");
    end PbBr;

    record PbBr2 "PbBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbBr2);
     annotation (preferredView = "info");
    end PbBr2;

    record PbBr3 "PbBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbBr3);
     annotation (preferredView = "info");
    end PbBr3;

    record PbBr4 "PbBr4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbBr4);
     annotation (preferredView = "info");
    end PbBr4;

    record PbCL "PbCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbCL);
     annotation (preferredView = "info");
    end PbCL;

    record PbCL2 "PbCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbCL2);
     annotation (preferredView = "info");
    end PbCL2;

    record PbCL3 "PbCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbCL3);
     annotation (preferredView = "info");
    end PbCL3;

    record PbCL4 "PbCL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbCL4);
     annotation (preferredView = "info");
    end PbCL4;

    record PbF "PbF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbF);
     annotation (preferredView = "info");
    end PbF;

    record PbF2 "PbF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbF2);
     annotation (preferredView = "info");
    end PbF2;

    record PbF3 "PbF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbF3);
     annotation (preferredView = "info");
    end PbF3;

    record PbF4 "PbF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbF4);
     annotation (preferredView = "info");
    end PbF4;

    record PbI "PbI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbI);
     annotation (preferredView = "info");
    end PbI;

    record PbI2 "PbI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbI2);
     annotation (preferredView = "info");
    end PbI2;

    record PbI3 "PbI3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbI3);
     annotation (preferredView = "info");
    end PbI3;

    record PbI4 "PbI4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbI4);
     annotation (preferredView = "info");
    end PbI4;

    record PbO "PbO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbO);
     annotation (preferredView = "info");
    end PbO;

    record PbO2 "PbO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbO2);
     annotation (preferredView = "info");
    end PbO2;

    record PbS "PbS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbS);
     annotation (preferredView = "info");
    end PbS;

    record PbS2 "PbS2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.PbS2);
     annotation (preferredView = "info");
    end PbS2;

    record Rb "Rb(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb);
     annotation (preferredView = "info");
    end Rb;

    record Rbplus "Rbplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rbplus,
      z=1);
     annotation (preferredView = "info");
    end Rbplus;

    record Rbminus "Rbminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rbminus,
      z=-1);
     annotation (preferredView = "info");
    end Rbminus;

    record RbBO2 "RbBO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbBO2);
     annotation (preferredView = "info");
    end RbBO2;

    record RbBr "RbBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbBr);
     annotation (preferredView = "info");
    end RbBr;

    record RbCL "RbCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbCL);
     annotation (preferredView = "info");
    end RbCL;

    record RbF "RbF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbF);
     annotation (preferredView = "info");
    end RbF;

    record RbH "RbH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbH);
     annotation (preferredView = "info");
    end RbH;

    record RbI "RbI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbI);
     annotation (preferredView = "info");
    end RbI;

    record RbK "RbK(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbK);
     annotation (preferredView = "info");
    end RbK;

    record RbLi "RbLi(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbLi);
     annotation (preferredView = "info");
    end RbLi;

    record RbNO2 "RbNO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbNO2);
     annotation (preferredView = "info");
    end RbNO2;

    record RbNO3 "RbNO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbNO3);
     annotation (preferredView = "info");
    end RbNO3;

    record RbNa "RbNa(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbNa);
     annotation (preferredView = "info");
    end RbNa;

    record RbO "RbO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbO);
     annotation (preferredView = "info");
    end RbO;

    record RbOH "RbOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.RbOH);
     annotation (preferredView = "info");
    end RbOH;

    record Rb2Br2 "Rb2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2Br2);
     annotation (preferredView = "info");
    end Rb2Br2;

    record Rb2CL2 "Rb2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2CL2);
     annotation (preferredView = "info");
    end Rb2CL2;

    record Rb2F2 "Rb2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2F2);
     annotation (preferredView = "info");
    end Rb2F2;

    record Rb2I2 "Rb2I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2I2);
     annotation (preferredView = "info");
    end Rb2I2;

    record Rb2O "Rb2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2O);
     annotation (preferredView = "info");
    end Rb2O;

    record Rb2O2 "Rb2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2O2);
     annotation (preferredView = "info");
    end Rb2O2;

    record Rb2O2H2 "Rb2O2H2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2O2H2);
     annotation (preferredView = "info");
    end Rb2O2H2;

    record Rb2SO4 "Rb2SO4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rb2SO4);
     annotation (preferredView = "info");
    end Rb2SO4;

    record Rn "Rn(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rn);
     annotation (preferredView = "info");
    end Rn;

    record Rnplus "Rnplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Rnplus,
      z=1);
     annotation (preferredView = "info");
    end Rnplus;

    record S "S(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S);
     annotation (preferredView = "info");
    end S;

    record Splus "Splus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Splus,
      z=1);
     annotation (preferredView = "info");
    end Splus;

    record Sminus "Sminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Sminus,
      z=-1);
     annotation (preferredView = "info");
    end Sminus;

    record SCL "SCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SCL);
     annotation (preferredView = "info");
    end SCL;

    record SCL2 "SCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SCL2);
     annotation (preferredView = "info");
    end SCL2;

    record SCL2plus "SCL2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SCL2plus,
      z=1);
     annotation (preferredView = "info");
    end SCL2plus;

    record SD "SD(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SD);
     annotation (preferredView = "info");
    end SD;

    record SF "SF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF);
     annotation (preferredView = "info");
    end SF;

    record SFplus "SFplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SFplus,
      z=1);
     annotation (preferredView = "info");
    end SFplus;

    record SFminus "SFminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SFminus,
      z=-1);
     annotation (preferredView = "info");
    end SFminus;

    record SF2 "SF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF2);
     annotation (preferredView = "info");
    end SF2;

    record SF2plus "SF2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF2plus,
      z=1);
     annotation (preferredView = "info");
    end SF2plus;

    record SF2minus "SF2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF2minus,
      z=-1);
     annotation (preferredView = "info");
    end SF2minus;

    record SF3 "SF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF3);
     annotation (preferredView = "info");
    end SF3;

    record SF3plus "SF3plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF3plus,
      z=1);
     annotation (preferredView = "info");
    end SF3plus;

    record SF3minus "SF3minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF3minus,
      z=-1);
     annotation (preferredView = "info");
    end SF3minus;

    record SF4 "SF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF4);
     annotation (preferredView = "info");
    end SF4;

    record SF4plus "SF4plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF4plus,
      z=1);
     annotation (preferredView = "info");
    end SF4plus;

    record SF4minus "SF4minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF4minus,
      z=-1);
     annotation (preferredView = "info");
    end SF4minus;

    record SF5 "SF5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF5);
     annotation (preferredView = "info");
    end SF5;

    record SF5plus "SF5plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF5plus,
      z=1);
     annotation (preferredView = "info");
    end SF5plus;

    record SF5minus "SF5minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF5minus,
      z=-1);
     annotation (preferredView = "info");
    end SF5minus;

    record SF6 "SF6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF6);
     annotation (preferredView = "info");
    end SF6;

    record SF6minus "SF6minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SF6minus,
      z=-1);
     annotation (preferredView = "info");
    end SF6minus;

    record SH "SH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SH);
     annotation (preferredView = "info");
    end SH;

    record SHminus "SHminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SHminus,
      z=-1);
     annotation (preferredView = "info");
    end SHminus;

    record SN "SN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SN);
     annotation (preferredView = "info");
    end SN;

    record SO "SO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SO);
     annotation (preferredView = "info");
    end SO;

    record SOminus "SOminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SOminus,
      z=-1);
     annotation (preferredView = "info");
    end SOminus;

    record SOF2 "SOF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SOF2);
     annotation (preferredView = "info");
    end SOF2;

    record SO2 "SO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2);
     annotation (preferredView = "info");
    end SO2;

    record SO2minus "SO2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2minus,
      z=-1);
     annotation (preferredView = "info");
    end SO2minus;

    record SO2CL2 "SO2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2CL2);
     annotation (preferredView = "info");
    end SO2CL2;

    record SO2FCL "SO2FCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2FCL);
     annotation (preferredView = "info");
    end SO2FCL;

    record SO2F2 "SO2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SO2F2);
     annotation (preferredView = "info");
    end SO2F2;

    record SO3 "SO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SO3);
     annotation (preferredView = "info");
    end SO3;

    record S2 "S2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S2);
     annotation (preferredView = "info");
    end S2;

    record S2minus "S2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S2minus,
      z=-1);
     annotation (preferredView = "info");
    end S2minus;

    record S2CL2 "S2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S2CL2);
     annotation (preferredView = "info");
    end S2CL2;

    record S2F2 "S2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S2F2);
     annotation (preferredView = "info");
    end S2F2;

    record S2O "S2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S2O);
     annotation (preferredView = "info");
    end S2O;

    record S3 "S3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S3);
     annotation (preferredView = "info");
    end S3;

    record S4 "S4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S4);
     annotation (preferredView = "info");
    end S4;

    record S5 "S5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S5);
     annotation (preferredView = "info");
    end S5;

    record S6 "S6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S6);
     annotation (preferredView = "info");
    end S6;

    record S7 "S7(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S7);
     annotation (preferredView = "info");
    end S7;

    record S8 "S8(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.S8);
     annotation (preferredView = "info");
    end S8;

    record Sc "Sc(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Sc);
     annotation (preferredView = "info");
    end Sc;

    record Scplus "Scplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Scplus,
      z=1);
     annotation (preferredView = "info");
    end Scplus;

    record Scminus "Scminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Scminus,
      z=-1);
     annotation (preferredView = "info");
    end Scminus;

    record ScO "ScO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ScO);
     annotation (preferredView = "info");
    end ScO;

    record ScOplus "ScOplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ScOplus,
      z=1);
     annotation (preferredView = "info");
    end ScOplus;

    record ScO2 "ScO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ScO2);
     annotation (preferredView = "info");
    end ScO2;

    record Sc2O "Sc2O(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Sc2O);
     annotation (preferredView = "info");
    end Sc2O;

    record Sc2O2 "Sc2O2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Sc2O2);
     annotation (preferredView = "info");
    end Sc2O2;

    record Si "Si(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Si);
     annotation (preferredView = "info");
    end Si;

    record Siplus "Siplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Siplus,
      z=1);
     annotation (preferredView = "info");
    end Siplus;

    record Siminus "Siminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Siminus,
      z=-1);
     annotation (preferredView = "info");
    end Siminus;

    record SiBr "SiBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiBr);
     annotation (preferredView = "info");
    end SiBr;

    record SiBr2 "SiBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiBr2);
     annotation (preferredView = "info");
    end SiBr2;

    record SiBr3 "SiBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiBr3);
     annotation (preferredView = "info");
    end SiBr3;

    record SiBr4 "SiBr4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiBr4);
     annotation (preferredView = "info");
    end SiBr4;

    record SiC "SiC(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiC);
     annotation (preferredView = "info");
    end SiC;

    record SiC2 "SiC2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiC2);
     annotation (preferredView = "info");
    end SiC2;

    record SiCL "SiCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiCL);
     annotation (preferredView = "info");
    end SiCL;

    record SiCL2 "SiCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiCL2);
     annotation (preferredView = "info");
    end SiCL2;

    record SiCL3 "SiCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiCL3);
     annotation (preferredView = "info");
    end SiCL3;

    record SiCL4 "SiCL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiCL4);
     annotation (preferredView = "info");
    end SiCL4;

    record SiF "SiF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiF);
     annotation (preferredView = "info");
    end SiF;

    record SiFCL "SiFCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiFCL);
     annotation (preferredView = "info");
    end SiFCL;

    record SiF2 "SiF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiF2);
     annotation (preferredView = "info");
    end SiF2;

    record SiF3 "SiF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiF3);
     annotation (preferredView = "info");
    end SiF3;

    record SiF4 "SiF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiF4);
     annotation (preferredView = "info");
    end SiF4;

    record SiH "SiH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH);
     annotation (preferredView = "info");
    end SiH;

    record SiHplus "SiHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHplus,
      z=1);
     annotation (preferredView = "info");
    end SiHplus;

    record SiHBr3 "SiHBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHBr3);
     annotation (preferredView = "info");
    end SiHBr3;

    record SiHCL "SiHCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHCL);
     annotation (preferredView = "info");
    end SiHCL;

    record SiHCL3 "SiHCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHCL3);
     annotation (preferredView = "info");
    end SiHCL3;

    record SiHF "SiHF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHF);
     annotation (preferredView = "info");
    end SiHF;

    record SiHF3 "SiHF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHF3);
     annotation (preferredView = "info");
    end SiHF3;

    record SiHI3 "SiHI3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiHI3);
     annotation (preferredView = "info");
    end SiHI3;

    record SiH2 "SiH2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2);
     annotation (preferredView = "info");
    end SiH2;

    record SiH2Br2 "SiH2Br2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2Br2);
     annotation (preferredView = "info");
    end SiH2Br2;

    record SiH2CL2 "SiH2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2CL2);
     annotation (preferredView = "info");
    end SiH2CL2;

    record SiH2F2 "SiH2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2F2);
     annotation (preferredView = "info");
    end SiH2F2;

    record SiH2I2 "SiH2I2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH2I2);
     annotation (preferredView = "info");
    end SiH2I2;

    record SiH3 "SiH3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3);
     annotation (preferredView = "info");
    end SiH3;

    record SiH3Br "SiH3Br(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3Br);
     annotation (preferredView = "info");
    end SiH3Br;

    record SiH3CL "SiH3CL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3CL);
     annotation (preferredView = "info");
    end SiH3CL;

    record SiH3F "SiH3F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3F);
     annotation (preferredView = "info");
    end SiH3F;

    record SiH3I "SiH3I(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH3I);
     annotation (preferredView = "info");
    end SiH3I;

    record SiH4 "SiH4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiH4);
     annotation (preferredView = "info");
    end SiH4;

    record SiI "SiI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiI);
     annotation (preferredView = "info");
    end SiI;

    record SiI2 "SiI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiI2);
     annotation (preferredView = "info");
    end SiI2;

    record SiN "SiN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiN);
     annotation (preferredView = "info");
    end SiN;

    record SiO "SiO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiO);
     annotation (preferredView = "info");
    end SiO;

    record SiO2 "SiO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiO2);
     annotation (preferredView = "info");
    end SiO2;

    record SiS "SiS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiS);
     annotation (preferredView = "info");
    end SiS;

    record SiS2 "SiS2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SiS2);
     annotation (preferredView = "info");
    end SiS2;

    record Si2 "Si2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Si2);
     annotation (preferredView = "info");
    end Si2;

    record Si2C "Si2C(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Si2C);
     annotation (preferredView = "info");
    end Si2C;

    record Si2F6 "Si2F6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Si2F6);
     annotation (preferredView = "info");
    end Si2F6;

    record Si2N "Si2N(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Si2N);
     annotation (preferredView = "info");
    end Si2N;

    record Si3 "Si3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Si3);
     annotation (preferredView = "info");
    end Si3;

    record Sn "Sn(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Sn);
     annotation (preferredView = "info");
    end Sn;

    record Snplus "Snplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Snplus,
      z=1);
     annotation (preferredView = "info");
    end Snplus;

    record Snminus "Snminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Snminus,
      z=-1);
     annotation (preferredView = "info");
    end Snminus;

    record SnBr "SnBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnBr);
     annotation (preferredView = "info");
    end SnBr;

    record SnBr2 "SnBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnBr2);
     annotation (preferredView = "info");
    end SnBr2;

    record SnBr3 "SnBr3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnBr3);
     annotation (preferredView = "info");
    end SnBr3;

    record SnBr4 "SnBr4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnBr4);
     annotation (preferredView = "info");
    end SnBr4;

    record SnCL "SnCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnCL);
     annotation (preferredView = "info");
    end SnCL;

    record SnCL2 "SnCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnCL2);
     annotation (preferredView = "info");
    end SnCL2;

    record SnCL3 "SnCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnCL3);
     annotation (preferredView = "info");
    end SnCL3;

    record SnCL4 "SnCL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnCL4);
     annotation (preferredView = "info");
    end SnCL4;

    record SnF "SnF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnF);
     annotation (preferredView = "info");
    end SnF;

    record SnF2 "SnF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnF2);
     annotation (preferredView = "info");
    end SnF2;

    record SnF3 "SnF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnF3);
     annotation (preferredView = "info");
    end SnF3;

    record SnF4 "SnF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnF4);
     annotation (preferredView = "info");
    end SnF4;

    record SnI "SnI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnI);
     annotation (preferredView = "info");
    end SnI;

    record SnI2 "SnI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnI2);
     annotation (preferredView = "info");
    end SnI2;

    record SnI3 "SnI3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnI3);
     annotation (preferredView = "info");
    end SnI3;

    record SnI4 "SnI4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnI4);
     annotation (preferredView = "info");
    end SnI4;

    record SnO "SnO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnO);
     annotation (preferredView = "info");
    end SnO;

    record SnO2 "SnO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnO2);
     annotation (preferredView = "info");
    end SnO2;

    record SnS "SnS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnS);
     annotation (preferredView = "info");
    end SnS;

    record SnS2 "SnS2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SnS2);
     annotation (preferredView = "info");
    end SnS2;

    record Sn2 "Sn2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Sn2);
     annotation (preferredView = "info");
    end Sn2;

    record Sr "Sr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Sr);
     annotation (preferredView = "info");
    end Sr;

    record Srplus "Srplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Srplus,
      z=1);
     annotation (preferredView = "info");
    end Srplus;

    record SrBr "SrBr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrBr);
     annotation (preferredView = "info");
    end SrBr;

    record SrBr2 "SrBr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrBr2);
     annotation (preferredView = "info");
    end SrBr2;

    record SrCL "SrCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrCL);
     annotation (preferredView = "info");
    end SrCL;

    record SrCLplus "SrCLplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrCLplus,
      z=1);
     annotation (preferredView = "info");
    end SrCLplus;

    record SrCL2 "SrCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrCL2);
     annotation (preferredView = "info");
    end SrCL2;

    record SrF "SrF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrF);
     annotation (preferredView = "info");
    end SrF;

    record SrFplus "SrFplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrFplus,
      z=1);
     annotation (preferredView = "info");
    end SrFplus;

    record SrF2 "SrF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrF2);
     annotation (preferredView = "info");
    end SrF2;

    record SrH "SrH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrH);
     annotation (preferredView = "info");
    end SrH;

    record SrI "SrI(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrI);
     annotation (preferredView = "info");
    end SrI;

    record SrI2 "SrI2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrI2);
     annotation (preferredView = "info");
    end SrI2;

    record SrO "SrO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrO);
     annotation (preferredView = "info");
    end SrO;

    record SrOH "SrOH(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrOH);
     annotation (preferredView = "info");
    end SrOH;

    record SrOHplus "SrOHplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrOHplus,
      z=1);
     annotation (preferredView = "info");
    end SrOHplus;

    record Sr_OH_2 "Sr_OH_2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Sr_OH_2);
     annotation (preferredView = "info");
    end Sr_OH_2;

    record SrS "SrS(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.SrS);
     annotation (preferredView = "info");
    end SrS;

    record Sr2 "Sr2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Sr2);
     annotation (preferredView = "info");
    end Sr2;

    record Ta "Ta(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ta);
     annotation (preferredView = "info");
    end Ta;

    record Taplus "Taplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Taplus,
      z=1);
     annotation (preferredView = "info");
    end Taplus;

    record Taminus "Taminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Taminus,
      z=-1);
     annotation (preferredView = "info");
    end Taminus;

    record TaCL5 "TaCL5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TaCL5);
     annotation (preferredView = "info");
    end TaCL5;

    record TaO "TaO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TaO);
     annotation (preferredView = "info");
    end TaO;

    record TaO2 "TaO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TaO2);
     annotation (preferredView = "info");
    end TaO2;

    record Ti "Ti(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Ti);
     annotation (preferredView = "info");
    end Ti;

    record Tiplus "Tiplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Tiplus,
      z=1);
     annotation (preferredView = "info");
    end Tiplus;

    record Timinus "Timinus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Timinus,
      z=-1);
     annotation (preferredView = "info");
    end Timinus;

    record TiCL "TiCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TiCL);
     annotation (preferredView = "info");
    end TiCL;

    record TiCL2 "TiCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TiCL2);
     annotation (preferredView = "info");
    end TiCL2;

    record TiCL3 "TiCL3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TiCL3);
     annotation (preferredView = "info");
    end TiCL3;

    record TiCL4 "TiCL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TiCL4);
     annotation (preferredView = "info");
    end TiCL4;

    record TiO "TiO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TiO);
     annotation (preferredView = "info");
    end TiO;

    record TiOplus "TiOplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TiOplus,
      z=1);
     annotation (preferredView = "info");
    end TiOplus;

    record TiOCL "TiOCL(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TiOCL);
     annotation (preferredView = "info");
    end TiOCL;

    record TiOCL2 "TiOCL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TiOCL2);
     annotation (preferredView = "info");
    end TiOCL2;

    record TiO2 "TiO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.TiO2);
     annotation (preferredView = "info");
    end TiO2;

    record U "U(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.U);
     annotation (preferredView = "info");
    end U;

    record UF "UF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF);
     annotation (preferredView = "info");
    end UF;

    record UFplus "UFplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UFplus,
      z=1);
     annotation (preferredView = "info");
    end UFplus;

    record UFminus "UFminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UFminus,
      z=-1);
     annotation (preferredView = "info");
    end UFminus;

    record UF2 "UF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF2);
     annotation (preferredView = "info");
    end UF2;

    record UF2plus "UF2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF2plus,
      z=1);
     annotation (preferredView = "info");
    end UF2plus;

    record UF2minus "UF2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF2minus,
      z=-1);
     annotation (preferredView = "info");
    end UF2minus;

    record UF3 "UF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF3);
     annotation (preferredView = "info");
    end UF3;

    record UF3plus "UF3plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF3plus,
      z=1);
     annotation (preferredView = "info");
    end UF3plus;

    record UF3minus "UF3minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF3minus,
      z=-1);
     annotation (preferredView = "info");
    end UF3minus;

    record UF4 "UF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF4);
     annotation (preferredView = "info");
    end UF4;

    record UF4plus "UF4plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF4plus,
      z=1);
     annotation (preferredView = "info");
    end UF4plus;

    record UF4minus "UF4minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF4minus,
      z=-1);
     annotation (preferredView = "info");
    end UF4minus;

    record UF5 "UF5(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF5);
     annotation (preferredView = "info");
    end UF5;

    record UF5plus "UF5plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF5plus,
      z=1);
     annotation (preferredView = "info");
    end UF5plus;

    record UF5minus "UF5minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF5minus,
      z=-1);
     annotation (preferredView = "info");
    end UF5minus;

    record UF6 "UF6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF6);
     annotation (preferredView = "info");
    end UF6;

    record UF6minus "UF6minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UF6minus,
      z=-1);
     annotation (preferredView = "info");
    end UF6minus;

    record UO "UO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UO);
     annotation (preferredView = "info");
    end UO;

    record UOplus "UOplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UOplus,
      z=1);
     annotation (preferredView = "info");
    end UOplus;

    record UOF "UOF(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UOF);
     annotation (preferredView = "info");
    end UOF;

    record UOF2 "UOF2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UOF2);
     annotation (preferredView = "info");
    end UOF2;

    record UOF3 "UOF3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UOF3);
     annotation (preferredView = "info");
    end UOF3;

    record UOF4 "UOF4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UOF4);
     annotation (preferredView = "info");
    end UOF4;

    record UO2 "UO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2);
     annotation (preferredView = "info");
    end UO2;

    record UO2plus "UO2plus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2plus,
      z=1);
     annotation (preferredView = "info");
    end UO2plus;

    record UO2minus "UO2minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2minus,
      z=-1);
     annotation (preferredView = "info");
    end UO2minus;

    record UO2F "UO2F(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2F);
     annotation (preferredView = "info");
    end UO2F;

    record UO2F2 "UO2F2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UO2F2);
     annotation (preferredView = "info");
    end UO2F2;

    record UO3 "UO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UO3);
     annotation (preferredView = "info");
    end UO3;

    record UO3minus "UO3minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.UO3minus,
      z=-1);
     annotation (preferredView = "info");
    end UO3minus;

    record V "V(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.V);
     annotation (preferredView = "info");
    end V;

    record Vplus "Vplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Vplus,
      z=1);
     annotation (preferredView = "info");
    end Vplus;

    record Vminus "Vminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Vminus,
      z=-1);
     annotation (preferredView = "info");
    end Vminus;

    record VCL4 "VCL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.VCL4);
     annotation (preferredView = "info");
    end VCL4;

    record VN "VN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.VN);
     annotation (preferredView = "info");
    end VN;

    record VO "VO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.VO);
     annotation (preferredView = "info");
    end VO;

    record VO2 "VO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.VO2);
     annotation (preferredView = "info");
    end VO2;

    record V4O10 "V4O10(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.V4O10);
     annotation (preferredView = "info");
    end V4O10;

    record W "W(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.W);
     annotation (preferredView = "info");
    end W;

    record Wplus "Wplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Wplus,
      z=1);
     annotation (preferredView = "info");
    end Wplus;

    record Wminus "Wminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Wminus,
      z=-1);
     annotation (preferredView = "info");
    end Wminus;

    record WCL6 "WCL6(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.WCL6);
     annotation (preferredView = "info");
    end WCL6;

    record WO "WO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.WO);
     annotation (preferredView = "info");
    end WO;

    record WOCL4 "WOCL4(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.WOCL4);
     annotation (preferredView = "info");
    end WOCL4;

    record WO2 "WO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.WO2);
     annotation (preferredView = "info");
    end WO2;

    record WO2CL2 "WO2CL2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.WO2CL2);
     annotation (preferredView = "info");
    end WO2CL2;

    record WO3 "WO3(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.WO3);
     annotation (preferredView = "info");
    end WO3;

    record WO3minus "WO3minus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.WO3minus,
      z=-1);
     annotation (preferredView = "info");
    end WO3minus;

    record Xe "Xe(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Xe);
     annotation (preferredView = "info");
    end Xe;

    record Xeplus "Xeplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Xeplus,
      z=1);
     annotation (preferredView = "info");
    end Xeplus;

    record Zn "Zn(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Zn);
     annotation (preferredView = "info");
    end Zn;

    record Znplus "Znplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Znplus,
      z=1);
     annotation (preferredView = "info");
    end Znplus;

    record Zr "Zr(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Zr);
     annotation (preferredView = "info");
    end Zr;

    record Zrplus "Zrplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Zrplus,
      z=1);
     annotation (preferredView = "info");
    end Zrplus;

    record Zrminus "Zrminus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.Zrminus,
      z=-1);
     annotation (preferredView = "info");
    end Zrminus;

    record ZrN "ZrN(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ZrN);
     annotation (preferredView = "info");
    end ZrN;

    record ZrO "ZrO(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ZrO);
     annotation (preferredView = "info");
    end ZrO;

    record ZrOplus "ZrOplus(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ZrOplus,
      z=1);
     annotation (preferredView = "info");
    end ZrOplus;

    record ZrO2 "ZrO2(g) MSL"
     extends Chemical.Obsolete.Interfaces.IdealGasMSL.SubstanceData
                                                          (
      data = Modelica.Media.IdealGases.Common.SingleGasesData.ZrO2);
     annotation (preferredView = "info");
    end ZrO2;
    end IdealGasesMSL;
      extends Modelica.Icons.Package;

    record Silver_solid "Ag(s)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.1078682,
        z=0,
        DfH=0,
        DfG=0,
        Cp=25.4,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"});
          annotation (preferredView = "info");
    end Silver_solid;

    record Silver_aqueous "Ag+(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.1078682,
        z=1,
        DfH=105900,
        DfG=77100,
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end Silver_aqueous;

    record SilverChloride_solid "AgCl(s)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.14332,
        z=0,
        DfH=-127030,
        DfG=-109720,
        Cp=50.8,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end SilverChloride_solid;

    record Calcium_aqueous "Ca++(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.0401,
        z=2,
        DfH=-542960,
        DfG=-542960 - 298.15*(33.67),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end Calcium_aqueous;

    record Chloride_aqueous "Cl-(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.03545,
        z=-1,
        DfH=-167460,
        DfG=-131170,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Chloride_aqueous;

    record CarbonMonoxide_gas "CO(g)"
     extends Chemical.Obsolete.Interfaces.IdealGas.SubstanceData
                                                       (
        MolarWeight=0.02801,
        DfH=-110500,
        DfG=-137300,
        Cp=29.13,
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.engineeringtoolbox.com/carbon-monoxide-d_975.html"});
          annotation (preferredView = "info");
    end CarbonMonoxide_gas;

    record CarbonMonoxide_aqueous "CO(aq*)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.02801,
        DfH=-276900,
        DfG=-110200,
        References={
            "Calculated from gas phase using Henry's coefficient from http://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Mask=10"});
      annotation (preferredView = "info");
    end CarbonMonoxide_aqueous;
              //  DfG = -8.314*298.15*log(0.00099/55.508)  +  -137300

    record CarbonDioxide_gas "CO2(g)"
     extends Chemical.Obsolete.Interfaces.IdealGas.SubstanceData
                                                       (
        MolarWeight=0.044,
        DfH=-393500,
        DfG=-394400,
        Cp=37.1,
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end CarbonDioxide_gas;

    record CarbonDioxide_aqueous "CO2(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        gamma=1.17385,
        MolarWeight=0.044,
        DfH=-412900,
        DfG=-386200,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end CarbonDioxide_aqueous;

    record Carbonate_aqueous "CO3--(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.06001,
        z=-2,
        DfH=-676300,
        DfG=-676300 - 298.15*(-497.065),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Carbonate_aqueous;

    record Electrone_solid "e-(s)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=5.4857990946e-7,
        z=-1,
        DfH=0,
        DfG=0,
        Cp=0,
        References={
            "http://physics.nist.gov/cgi-bin/cuu/Value?mme, To solve standard electo-chemical cell potentials"});
          annotation (preferredView = "info");
    end Electrone_solid;

    record Iron2_aqueous "Fe++(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.05585,
        z=2,
        DfH=-87860,
        DfG=-87860 - 298.15*(-9.93),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end Iron2_aqueous;

    record Iron3_aqueous "Fe+++(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.05585,
        z=3,
        DfH=-47700,
        DfG=-47700 - 298.15*(-124.77),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end Iron3_aqueous;

    record Glucose_solid "Glu(s)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.1806,
        DfH=-1274500,
        DfG=-1274500 - 298.15*(-1220.66),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end Glucose_solid;

    record Hydrogen_gas "H2(g)"
     extends Chemical.Obsolete.Interfaces.IdealGas.SubstanceData
                                                       (
        MolarWeight=0.00201588,
        z=0,
        DfH=0,
        DfG=0,
        Cp=28.8,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Hydrogen_gas;

    record CarbonicAcid_aqueous "H2CO3(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.062027,
        DfH=-699700,
        DfG=-699700 - 298.15*(-256.582),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end CarbonicAcid_aqueous;

    record Water_gas "H2O(g)"
     extends Chemical.Obsolete.Interfaces.IdealGas.SubstanceData
                                                       (
        MolarWeight=0.018015,
        DfH=-241830,
        DfG=-228590,
        Cp=33.6,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Water_gas;

    record Water_liquid_without_selfClustering "H2O(l) without self-clustering"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.018015,
        DfH=-285840,
        DfG=-237190,
        Cp=75.3,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});

        annotation (preferredView = "info", Documentation(info="<html>
<p><br><span style=\"font-family: Courier New;\">&nbsp;&nbsp;&nbsp;&nbsp;</span></p>
</html>"));
    end Water_liquid_without_selfClustering;
    //   Cv=74.539,
    // Enthalpy as in H2O(l) = with assumption that hydrogen bonds do not have significant enthaplies

    record Water_liquid "H2O(l) with self-clustering"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.018015,
        DfH=-285830,
        DfG=-227230,
        Cp=75.3,
        SelfClustering = true,
        SelfClustering_dH = -81.6348,
        SelfClustering_dS = 32.845554,
          References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
                                      //-77.95928,

    /*  SelfClustering_dH = -81.6348,
    SelfClustering_dS = 32.8344,
*/

      // S=(0 + Modelica.Constants.R*(273.15+25)*log(55.345/0.95-1))/(273.15+25),
      // SelfClustering_dS = (SelfClustering_dH + Modelica.Constants.R*(273.15+25)*log((55.345-1)/1))/(273.15+25),
      annotation (preferredView = "info", Documentation(info="<html>
<p><span style=\"font-family: Courier New;\">Even the tabulated formation Gibbs energy is DfG=-237190 there is another values because of water self-clustering. </span></p>
<p><br><span style=\"font-family: Courier New;\">New reported values for free water molecule in solution is calculated from water dissociation reaction.</span></p>
<p><span style=\"font-family: Courier New;\">&nbsp;&nbsp;&nbsp;&nbsp;</span></p>
</html>"));
    end Water_liquid;

    record Water_IceIh "H2O(s) - Ice I h"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.018015,
        DfH=-292639,
        DfG=-236590,
        Cp=37.77,
        References={"http://www1.lsbu.ac.uk/water/water_properties.html#pot"});
      annotation (preferredView = "info");
    end Water_IceIh;

    record DihydrogenPhosphate_aqueous "H2PO4-(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.095,
        z=-1,
        DfH=-1302480,
        DfG=-1302480 - 298.15*(-561.395),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end DihydrogenPhosphate_aqueous;

    record Hydronium_aqueous "H3O+(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.019022,
        z=1,
        DfH=-285840,
        DfG=-285840 - 298.15*(-163.17),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Hydronium_aqueous;

    record PhosphoricAcid_aqueous "H3PO4(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.095,
        DfH=-1288000,
        DfG=-1288000 - 298.15*(-496.4),
        References={
            "https://en.wikipedia.org/wiki/Phosphoric_acid, https://www.researchgate.net/publication/6600409_Standard_thermodynamic_properties_of_H3PO4%28aq%29_over_a_wide_range_of_temperatures_and_pressures"});
          annotation (preferredView = "info");
    end PhosphoricAcid_aqueous;

    record Proton_aqueous "H+(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.001007,
        z=1,
        DfH=0,
        DfG=0,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Proton_aqueous;
                       // as hypothetical HA <-> H+ + A- simplification of H2O + HA <-> H3O+ + A-";

    record Bicarbonate_aqueous "HCO3-(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.06102,
        z=-1,
        DfH=-691100,
        DfG=-691100 - 298.15*(-348.82),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Bicarbonate_aqueous;

    record Bicarbonate_blood "HCO3-(blood)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.06102,
        z=-1,
        DfH=-691100,
        DfG=-691100 - 298.15*(-348.82),
        gamma=0.79,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Bicarbonate_blood;

    record HydrogenPhosphate_aqueous "HPO4--(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.095,
        z=-2,
        DfH=-1298700,
        DfG=-1298700 - 298.15*(-686.232),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end HydrogenPhosphate_aqueous;

    record HydrogenSulfate_aqueous "HSO4-(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.097,
        z=-1,
        DfH=-885750,
        DfG=-752870,
        density=1800,
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end HydrogenSulfate_aqueous;

    record Potassium_aqueous "K+(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.0391,
        z=1,
        DfH=-251200,
        DfG=-251200 - 298.15*(103.97),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Potassium_aqueous;

    record Magnesium_aqueous "Mg++(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.0243,
        z=2,
        DfH=-461960,
        DfG=-461960 - 298.15*(-19.99),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.vias.org/genchem/standard_enthalpies_table.html"});
          annotation (preferredView = "info");
    end Magnesium_aqueous;

    record Sodium_aqueous "Na+(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.02299,
        z=1,
        DfH=-239660,
        DfG=-239660 - 298.15*(74.49),
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Sodium_aqueous;

    record Amonium_aqueous "NH4+(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.01804,
        z=1,
        DfH=-132800,
        DfG=-132800 - 298.15*(-178.77),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end Amonium_aqueous;

    record Oxygen_gas "O2(g)"
     extends Chemical.Obsolete.Interfaces.IdealGas.SubstanceData
                                                       (
        MolarWeight=0.032,
        DfH=0,
        DfG=0,
        Cp=29.4,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Oxygen_gas;

    record Oxygen_gas_Shomate_298_6000 "O2(g) Shomate 298K–6000K"
     extends Chemical.Obsolete.Interfaces.IdealGasShomate.SubstanceData
                                                              (
        MolarWeight=0.032,
        DfH=0,
        DfG=0,
        Cp=29.4,
        B=6.137261,
        C=-1.186521,
        D=0.09578,
        E=-0.219663,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html, http://old.vscht.cz/fch/cz/pomucky/fchab/C.html"});
      annotation (preferredView = "info");
    end Oxygen_gas_Shomate_298_6000;

    record Oxygen_gas_Shomate_200_5000 "O2(g) Shomate 200K–5000K"
     extends Chemical.Obsolete.Interfaces.IdealGasShomate.SubstanceData
                                                              (
        MolarWeight=0.032,
        DfH=0,
        DfG=0,
        Cp=29.4,
        B=-21.55543,
        C=2.456517,
        D=-0.16151,
        E=0.175056,
        X=44.837013,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html, http://old.vscht.cz/fch/cz/pomucky/fchab/C.html"});
      annotation (preferredView = "info");
    end Oxygen_gas_Shomate_200_5000;
              //A=8.99044,

    record Oxygen_aqueous "O2(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.032,
        DfH=-11700,
        DfG=16320,
        References={
            "http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.pdf, https://books.google.cz/books?id=dr-VBAAAQBAJ&pg=PA156&lpg=PA156&dq=Gibbs+energy+of+formation++%22O2(aq)%22&source=bl&ots=09N5CxY7OD&sig=hbsTXQvX59vXBqHUjFVVIZQpHCA&hl=cs&sa=X&ei=sDQtVaeUMMaRsAHpzYHgAg&redir_esc=y#v=onepage&q=Gibbs%20energy%20of%20formation%20%20%22O2(aq)%22&f=false"});
          annotation (preferredView = "info");
    end Oxygen_aqueous;

    record Hydroxide_aqueous "OH-(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.017006,
        z=-1,
        DfH=-229940,
        DfG=-157300,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Hydroxide_aqueous;

    record Lead_solid "Pb(s)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.2072,
        z=0,
        DfH=0,
        DfG=0,
        Cp=26.4,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"});
          annotation (preferredView = "info");
    end Lead_solid;

    record LeadDioxide_solid "PbO2(s)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.2391988,
        z=0,
        DfH=-276600,
        DfG=-219000,
        Cp=64.6,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"});
          annotation (preferredView = "info");
    end LeadDioxide_solid;

    record LeadSulfate_solid "PbSO4(s)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.30326,
        z=0,
        DfH=-918400,
        DfG=-811200,
        Cp=103.2,
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"});
          annotation (preferredView = "info");
    end LeadSulfate_solid;

    record Phosphate_aqueous "PO4---(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.095,
        z=-3,
        DfH=-1284070,
        DfG=-1284070 - 298.15*(-866.946),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end Phosphate_aqueous;

    record Sulphates_aqueous "SO4--(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.09607,
        z=-2,
        DfH=-907500,
        DfG=-907500 - 298.15*(-555.123),
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"});
          annotation (preferredView = "info");
    end Sulphates_aqueous;

    record Ethanol_liquid "Ethanol C2H5OH(l)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.04607,
        z=0,
        DfH=-276980,
        DfG=-174180,
        Cp=112.4,
        density=789,
        References={
            "http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, https://en.wikipedia.org/wiki/Ethanol_(data_page)"});
          annotation (preferredView = "info");
    end Ethanol_liquid;
      //Some organic molecules: https://www.e-education.psu.edu/drupal6/files/be497b/pdf/Bioenergetics_AppA.pdf
    //Other source: http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf

    record Urea_aqueous "Urea(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.06006,
        z=0,
        DfH=-333189,
        DfG=-197150,
        References={"https://en.wikipedia.org/wiki/Urea"});
          annotation (preferredView = "info");
    end Urea_aqueous;

    record Globulins_aqueous "Glb(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=66.5,
        z=-4,
        DfH=0,
        DfG=0,
        References={"https://en.wikipedia.org/wiki/Human_serum_albumin"});
          annotation (preferredView = "info");
    end Globulins_aqueous;

    record Albumin_aqueous "Alb(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=66.5,
        z=-11.4,
        DfH=0,
        DfG=0,
        References={"https://en.wikipedia.org/wiki/Human_serum_albumin"});
          annotation (preferredView = "info");
    end Albumin_aqueous;

    record ADP3_aqueous "ADP3(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.427201,
        z=-3,
        DfH=0,
        DfG=0,
        References={"relative - designed only for ATP hydrolysis example"});
      annotation (preferredView = "info");
    end ADP3_aqueous;

    record ATP4_aqueous "ATP^4-(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.427201,
        z=-4,
        DfH = -1.0263e+6,
        DfG = -882161,
        References={"relative - designed only for ATP hydrolysis example"});
        // dle reakce: ATP + H2O <-> ADP + H2PO4-  (G=-30.5 kJ/mol. H=-20 kJ/mol)
      annotation (preferredView = "info");
    end ATP4_aqueous;

    record ATP3_aqueous "ATP^3-(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.427201,
        z=-4,
        DfH = -1.0263e+6,
        DfG = -919245,
        References={"relative - designed only for ATP hydrolysis example"});
        // dle reakce: ATP^4- + H+ <-> ATP^3-  (pKa=6.5, H=0 kJ/mol)
      annotation (preferredView = "info");
    end ATP3_aqueous;

    model OxygenGasOnTemperature
      Real cp1,cp2;
      Real H1,H2;
      Real S1,S2;
      Real T;
    equation
      T=200+time;
      cp1 =Chemical.Obsolete.Interfaces.IdealGasShomate.molarHeatCapacityCp(Oxygen_gas_Shomate_298_6000(), T);
      cp2 =Chemical.Obsolete.Interfaces.IdealGasShomate.molarHeatCapacityCp(Oxygen_gas_Shomate_200_5000(), T);
      H1 =Chemical.Obsolete.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(Oxygen_gas_Shomate_298_6000(), T);
      H2 =Chemical.Obsolete.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(Oxygen_gas_Shomate_200_5000(), T);
      S1 =Chemical.Obsolete.Interfaces.IdealGasShomate.molarEntropyPure(Oxygen_gas_Shomate_298_6000(), T);
      S2 =Chemical.Obsolete.Interfaces.IdealGasShomate.molarEntropyPure(Oxygen_gas_Shomate_200_5000(), T);
    end OxygenGasOnTemperature;

    record Nitrogen_gas "N2(g)"
       extends Chemical.Obsolete.Interfaces.IdealGas.SubstanceData
                                                         (
          MolarWeight=0.0280134,
          DfH=0,
          DfG=0,
          Cp=29.1,
          References={
              "http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Type=JANAFG&Plot=on"});
        annotation (preferredView = "info");
    end Nitrogen_gas;

    record Methan_gas "CH4(g)"
     extends Chemical.Obsolete.Interfaces.IdealGas.SubstanceData
                                                       (
        MolarWeight=0.01604246,
        z=0,
        DfH = -74848,
        DfG = -50794,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"});

      annotation (preferredView = "info");
    end Methan_gas;

    record Methan_aqueous "CH4(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.01604246,
        z=0,
        DfH = -88151,
        DfG = -34504,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=10#Solubility"});

      annotation (preferredView = "info");
    end Methan_aqueous;

    record AceticAcid_gas "CH3COOH(g)"
     extends Chemical.Obsolete.Interfaces.IdealGas.SubstanceData
                                                       (
        MolarWeight=0.060052,
        z=0,
        DfH = -436071,
        DfG = -378978,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C64197&Mask=10#Solubility"});

      annotation (preferredView = "info");
    end AceticAcid_gas;

    record AceticAcid_aqueous "CH3COOH(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.060052,
        z=0,
        DfH = -488453,
        DfG = -399600,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"});

      annotation (preferredView = "info");
    end AceticAcid_aqueous;

    record Acetate_aqueous "CH3COO-(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.059052,
        z=-1,
        DfH = -488871,
        DfG = -372500,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"});
      annotation (preferredView = "info");
    end Acetate_aqueous;

    record Hydrogen_aqueous "H2(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.00201588,
        z=0,
        DfH=-4157,
        DfG=17740,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=10#Solubility"});
      annotation (preferredView = "info");
    end Hydrogen_aqueous;

    record Ethanol_gas "C2H5OH(g)"
     extends Chemical.Obsolete.Interfaces.IdealGas.SubstanceData
                                                       (
        MolarWeight=0.04607,
        z=0,
        DfH = -235400,
        DfG = -168600,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"});

      annotation (preferredView = "info");
    end Ethanol_gas;

    record Ethanol_aqueous "C2H5OH(aq)"
     extends Chemical.Obsolete.Interfaces.Incompressible.SubstanceData
                                                             (
        MolarWeight=0.04607,
        z=0,
        DfH=-290276,
        DfG=-181607,
        References={
            "http://www.vias.org/genchem/standard_enthalpies_table.html, https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Units=SI&Mask=10#Solubility"});
      annotation (preferredView = "info");
    end Ethanol_aqueous;
  end Substances;

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

      Modelica.Blocks.Interfaces.RealInput molalityInput(start=Molality,final unit="mol/kg")=n/KG
        if useMolalityInput
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

      Modelica.Blocks.Interfaces.RealInput molarConcentrationInput(start=MolarConcentration,final unit="mol/m3", displayUnit="mol/l")=n/L
        if useMolarityInput
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
        start=MoleFraction)=x
        if useMoleFractionInput annotation (HideResult=true, Placement(transformation(
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

      Modelica.Blocks.Interfaces.RealInput uInput(final unit="J/mol")=port_a.u
        if usePotentialInput annotation (HideResult=true, Placement(transformation(
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
          start=BufferValue)=bufferValue
          if useBufferValueInput annotation (HideResult=true, Placement(transformation(
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
        Chemical.Obsolete.Interfaces.Incompressible                         constrainedby Chemical.Obsolete.Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Obsolete.Interfaces.Incompressible
                                                "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Obsolete.Interfaces.IdealGas
                                                "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Obsolete.Interfaces.IdealGasMSL
                                                "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Obsolete.Interfaces.IdealGasShomate
                                                "Ideal Gas using Shomate model")));

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
            Chemical.Obsolete.Interfaces.Incompressible
                                                "Incompressible"),
          choice(redeclare package stateOfMatter =
            Chemical.Obsolete.Interfaces.IdealGas
                                                "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Obsolete.Interfaces.IdealGasMSL
                                                "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Obsolete.Interfaces.IdealGasShomate
                                                "Ideal Gas using Shomate model")));

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
        input Modelica.Units.SI.Mass massH2O=1 "Mass of H2O";
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
        input Modelica.Units.SI.Mass massH2O=1 "Mass of H2O in solution";
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

          amountOfBaseMolecules:=massH2O/substanceData.MolarWeight;
          x:=((2*SelfClustering_K+nSolution/amountOfBaseMolecules) -
           sqrt((4*SelfClustering_K*nSolution/amountOfBaseMolecules)+
           (nSolution/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

          specificAmountOfFreeBaseMolecule := (x*nSolution)/massH2O;

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
          Chemical.Obsolete.Interfaces.StateOfMatter constrainedby StateOfMatter
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
                Chemical.Obsolete.Interfaces.Incompressible
                                                    "Incompressible"),
              choice(redeclare package stateOfMatter =
                Chemical.Obsolete.Interfaces.IdealGas
                                                    "Ideal Gas"),
              choice(redeclare package stateOfMatter =
                Chemical.Obsolete.Interfaces.IdealGasMSL
                                                    "Ideal Gas from MSL"),
              choice(redeclare package stateOfMatter =
                Chemical.Obsolete.Interfaces.IdealGasShomate
                                                    "Ideal Gas using Shomate model")));

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

  package Examples "Examples that demonstrate usage of chemical library"
  extends Modelica.Icons.ExamplesPackage;

    model SimpleReaction
      "The simple chemical reaction A<->B with equilibrium B/A = 2"
       extends Modelica.Icons.Example;

      constant Real K = 2 "Dissociation constant of the reaction";

      constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Chemical.Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Obsolete.Components.Substance A(
        substanceData(MolarWeight=1),
        use_mass_start=false,
        amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-52,-8},{-32,12}})));

      Chemical.Obsolete.Components.Reaction reaction(nS=1, nP=1) annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
      Chemical.Obsolete.Components.Substance B(
        substanceData(DfG=-R*T_25degC*log(K), MolarWeight=1),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{62,-8},{42,12}})));

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
</html>",   info="<html>
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

      Chemical.Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Obsolete.Components.Substance A(use_mass_start=false, amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
      Chemical.Obsolete.Components.Reaction reaction(nS=2, nP=1) annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Chemical.Obsolete.Components.Substance B(use_mass_start=false, amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Chemical.Obsolete.Components.Substance C(
        substanceData(DfG=-R*T_25degC*log(Kx)),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{68,-8},{48,12}})));

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
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
        experiment(StopTime=0.0001, __Dymola_Algorithm="Dassl"));
    end SimpleReaction2;

    model HeatingOfWater "Heating of 1 kg water"
      extends Modelica.Icons.Example;

      Chemical.Obsolete.Components.Solution solution(useMechanicPorts=true, useThermalPort=true)
        annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
        annotation (Placement(transformation(extent={{-86,-72},{-66,-52}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed1
        annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
      Obsolete.Components.Substance liquidWater(mass_start=1, substanceData=Chemical.Obsolete.Substances.Water_liquid())
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
</html>",   info="<html>
<p>Heating of solution by one degree, using standard HeatPort from Modelica Standard Library.</p>
<p>Observe Solution.T (or H2O.Solution.T) for temperature change.</p>
</html>"));
    end HeatingOfWater;

    model HeatingOfAlcohol "Heating of 50% ethanol"
      extends Modelica.Icons.Example;

      Chemical.Obsolete.Components.Solution solution(useMechanicPorts=true, useThermalPort=true)
        annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
        annotation (Placement(transformation(extent={{-86,-76},{-66,-56}})));
      Chemical.Obsolete.Components.Substance Ethanol(
        redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
        substanceData=Chemical.Obsolete.Substances.Ethanol_liquid(),
        mass_start=(55.508/2)*0.04607) annotation (Placement(transformation(extent={{18,-8},{38,12}})));

      Modelica.Mechanics.Translational.Components.Fixed fixed1
        annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
      Obsolete.Components.Substance liquidWater(mass_start=1/2, substanceData=Chemical.Obsolete.Substances.Water_liquid())
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
</html>",   info="<html>
<p>Heating of solution of water and ethanol, using standard HeatPort from Modelica Standard Library.</p>
<p>Observe Solution.T (or H2O.Solution.T or Ethanol.Solution.T) for temperature change. Note, that we can heat all substances in solution at once and the results would differ from HeatingOfWater.</p>
</html>"));
    end HeatingOfAlcohol;

    model ExothermicReaction
      "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"

       extends Modelica.Icons.Example;

      parameter Modelica.Units.SI.MolarEnergy ReactionEnthalpy=-55000;

      Chemical.Obsolete.Components.Solution thermal_isolated_solution(useMechanicPorts=true, ConstantTemperature=false)
        annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));
      Chemical.Obsolete.Components.Substance A(use_mass_start=false, amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
      Chemical.Obsolete.Components.Reaction reaction(nS=1, nP=1) annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
      Chemical.Obsolete.Components.Substance B(
        substanceData(DfH=ReactionEnthalpy),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{40,-60},{20,-40}})));

      Chemical.Obsolete.Components.Solution solution_at_constant_temperature(useMechanicPorts=true, useThermalPort=true)
        annotation (Placement(transformation(extent={{-100,0},{98,94}})));
      Chemical.Obsolete.Components.Substance A1(use_mass_start=false, amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      Chemical.Obsolete.Components.Reaction reaction1(nS=1, nP=1) annotation (Placement(transformation(extent={{-8,40},{12,60}})));
      Chemical.Obsolete.Components.Substance B1(
        substanceData(DfH=ReactionEnthalpy),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{40,40},{20,60}})));

      //  Modelica.SIunits.HeatFlowRate q
      //    "Heat flow to environment to reach constant temperature";
      Modelica.Units.SI.Temperature t
        "Temperature if the solution is ideally thermal isolated from environment";
      Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{20,4},{40,24}})));
      Obsolete.Components.Substance H2O1(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
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
</html>",   info="<html>
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
      Chemical.Obsolete.Components.Solution idealGas(
        SurfaceArea=A,
        useMechanicPorts=true,
        useThermalPort=true,
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGas) annotation (Placement(transformation(extent={{-108,-50},{-8,50}})));
                       // AmbientPressure=p)
      //  volume_start=V,
      Chemical.Obsolete.Components.Substance H2_gas(
        redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.IdealGas,
        substanceData=Chemical.Obsolete.Substances.Hydrogen_gas(),
        use_mass_start=false,
        amountOfSubstance_start=0.026) annotation (Placement(transformation(extent={{-98,-26},{-78,-6}})));
      Chemical.Obsolete.Components.Substance O2_gas(
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=0.013) annotation (Placement(transformation(extent={{-98,10},{-78,30}})));
      Chemical.Obsolete.Components.Substance H2O_gas(
        substanceData=Chemical.Obsolete.Substances.Water_gas(),
        redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{-14,-8},{-34,12}})));
      Chemical.Obsolete.Components.Reaction reaction(
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
      Obsolete.Components.Solution idealGas1(
        SurfaceArea=A,
        useMechanicPorts=true,
        useThermalPort=true,
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGas) annotation (Placement(transformation(extent={{18,-52},{118,48}})));
      Obsolete.Components.Substance H2_gas1(
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGasMSL,
        substanceData(data=Modelica.Media.IdealGases.Common.SingleGasesData.H2),
        use_mass_start=false,
        amountOfSubstance_start=0.026) annotation (Placement(transformation(extent={{28,-28},{48,-8}})));
      Obsolete.Components.Substance O2_gas1(
        substanceData(data=Modelica.Media.IdealGases.Common.SingleGasesData.O2),
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGasMSL,
        use_mass_start=false,
        amountOfSubstance_start=0.013) annotation (Placement(transformation(extent={{28,8},{48,28}})));
      Obsolete.Components.Substance H2O_gas1(
        substanceData(data=Modelica.Media.IdealGases.Common.SingleGasesData.H2O),
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGasMSL,
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{112,-10},{92,10}})));

      Obsolete.Components.Reaction reaction1(
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

      Chemical.Obsolete.Components.Solution liquid(temperature_start=T_start, useThermalPort=true)
        annotation (Placement(transformation(extent={{-98,-98},{-6,-8}})));

      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
      Chemical.Obsolete.Components.Solution gas(
        temperature_start=T_start,
        useThermalPort=true,
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGas) annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                    /*volume_start(
        displayUnit="l") = 0.001, */
      Obsolete.Components.Substance H2O_gaseuous(
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGas,
        substanceData=Chemical.Obsolete.Substances.Water_gas(),
        mass_start=0.000106537) annotation (Placement(transformation(extent={{28,50},{8,70}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                           fixedTemperature
                  annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          origin={84,8})));
      Modelica.Blocks.Sources.ContinuousClock clock(offset=1*T_start)
        annotation (Placement(transformation(extent={{62,36},{82,56}})));
      Chemical.Obsolete.Components.Substance otherSubstances(
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGas,
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{2,28},{22,48}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(
        G=1e6) annotation (Placement(transformation(extent={{48,-8},{68,12}})));
      Obsolete.Components.Substance liquidWater(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{-28,-62},{-48,-42}})));
      inner Modelica.Fluid.System system(p_ambient=100000, T_ambient=298.15)
        annotation (Placement(transformation(extent={{54,-48},{74,-28}})));
      Obsolete.Components.GasSolubility gasSolubility annotation (Placement(transformation(extent={{-92,16},{-72,36}})));
      Obsolete.Sensors.PartialPressureSensor pH2O(redeclare package stateOfMatter = Obsolete.Interfaces.IdealGas, substanceData=
            Chemical.Obsolete.Substances.Water_gas()) annotation (Placement(transformation(extent={{-26,76},{-6,96}})));
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
</html>", revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end WaterVaporization;

    model WaterSublimation "Sublimation of water"
       extends Modelica.Icons.Example;

      parameter Modelica.Units.SI.Temperature T_start=273.15 - 50
        "Initial temperature";

      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
      Chemical.Obsolete.Components.Solution gas(
        temperature_start=T_start,
        useThermalPort=true,
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGas,
        BasePressure=600) annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                    /*volume_start(
        displayUnit="l") = 0.001, */
      Chemical.Obsolete.Components.Substance H2O_gaseuous(
        redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.IdealGas,
        substanceData=Chemical.Obsolete.Substances.Water_gas(),
        use_mass_start=false,
        amountOfSubstance_start=0.001) annotation (Placement(transformation(extent={{24,56},{4,76}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                           fixedTemperature
                  annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          origin={84,8})));
      Modelica.Blocks.Sources.ContinuousClock clock(offset=1*T_start)
        annotation (Placement(transformation(extent={{62,36},{82,56}})));
      Chemical.Obsolete.Components.Substance otherSubstances(
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        redeclare package stateOfMatter = Obsolete.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{-4,36},{16,56}})));
      Chemical.Obsolete.Components.Solution solid(
        temperature_start=T_start,
        BasePressure=600,
        useThermalPort=true) annotation (Placement(transformation(extent={{8,-98},{100,-8}})));
      Chemical.Obsolete.Components.Substance H2O_solid(
        substanceData=Chemical.Obsolete.Substances.Water_IceIh(),
        use_mass_start=false,
        amountOfSubstance_start=55.508) "Solid water" annotation (Placement(transformation(extent={{70,-62},{50,-42}})));
      Chemical.Obsolete.Components.GasSolubility gasSolubility1(KC=10) annotation (Placement(transformation(extent={{-76,18},{-56,38}})));
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
</html>", revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end WaterSublimation;

    model GasSolubility_NIST "Dissolution of gases in liquids"
       extends Modelica.Icons.Example;

      Chemical.Obsolete.Components.Solution water_solution_25degC(temperature_start=298.15) annotation (Placement(transformation(extent={{-160,-78},{-68,12}})));
                                          //(amountOfSolution_start=52.3)
      Chemical.Obsolete.Components.Solution water_solution_37degC(temperature_start=310.15) annotation (Placement(transformation(extent={{-52,-80},{42,12}})));
                                       //(amountOfSolution_start=39.7)
      Chemical.Obsolete.Components.GasSolubility CO2_dissolutionP annotation (Placement(transformation(extent={{-138,42},{-118,62}})));
      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
      Chemical.Obsolete.Components.Substance CO2_25(
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
        use_mass_start=false) "Free dissolved CO2 in water at 25 degC" annotation (Placement(transformation(extent={{-150,-26},{-130,-6}})));
      Chemical.Obsolete.Components.GasSolubility O2_dissolutionP annotation (Placement(transformation(extent={{-100,42},{-80,62}})));

      Chemical.Obsolete.Sources.ExternalIdealGasSubstance O2_g_25(
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        TotalPressure=system.p_ambient,
        PartialPressure(displayUnit="mmHg") = 12665.626804425,
        Temperature=298.15) annotation (Placement(transformation(extent={{-114,74},{-94,94}})));
      Chemical.Obsolete.Components.Substance O2_25(
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        substanceData=Chemical.Obsolete.Substances.Oxygen_aqueous(),
        use_mass_start=false) "Free dissolved O2 in water at 25 degC" annotation (Placement(transformation(extent={{-114,-26},{-94,-6}})));
      Chemical.Obsolete.Components.GasSolubility CO2_dissolutionE annotation (Placement(transformation(extent={{-26,42},{-6,62}})));

      Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO2_g_25(
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
        TotalPressure=system.p_ambient,
        PartialPressure(displayUnit="mmHg") = 5332.8954966,
        Temperature=298.15) annotation (Placement(transformation(extent={{-154,74},{-134,94}})));

      Chemical.Obsolete.Components.Substance CO2_37(
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
        use_mass_start=false) "Free dissolved CO2 in water at 37degC" annotation (Placement(transformation(extent={{-42,-34},{-22,-14}})));

      Chemical.Obsolete.Components.GasSolubility O2_dissolutionE_NIST annotation (Placement(transformation(extent={{18,42},{38,62}})));
      Chemical.Obsolete.Components.Substance O2_37(
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        substanceData=Chemical.Obsolete.Substances.Oxygen_aqueous(),
        use_mass_start=false) "Free dissolved O2 in water at 37degC" annotation (Placement(transformation(extent={{-2,-34},{18,-14}})));
      Obsolete.Components.Substance water_25(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{-100,-68},{-80,-48}})));
      Obsolete.Components.Substance water_37(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{8,-70},{28,-50}})));
      Obsolete.Sources.ExternalIdealGasSubstance CO2_g_37(
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
        TotalPressure=system.p_ambient,
        PartialPressure(displayUnit="mmHg") = 5332.8954966,
        Temperature=310.15) annotation (Placement(transformation(extent={{-44,68},{-24,88}})));
      Obsolete.Sources.ExternalIdealGasSubstance O2_g_37(
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        TotalPressure=system.p_ambient,
        PartialPressure(displayUnit="mmHg") = 12665.626804425,
        Temperature=310.15) annotation (Placement(transformation(extent={{-6,68},{14,88}})));
      Obsolete.Components.Solution water_solution_37degC1(temperature_start=273.15) annotation (Placement(transformation(extent={{66,-80},{160,12}})));
      Obsolete.Components.GasSolubility CO2_dissolutionE1 annotation (Placement(transformation(extent={{92,42},{112,62}})));
      Obsolete.Components.Substance CO2_0(
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
        use_mass_start=false) "Free dissolved CO2 in water at 0degC" annotation (Placement(transformation(extent={{76,-34},{96,-14}})));
      Obsolete.Components.GasSolubility O2_dissolutionE_NIST1 annotation (Placement(transformation(extent={{136,42},{156,62}})));
      Obsolete.Components.Substance O2_0(
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        substanceData=Chemical.Obsolete.Substances.Oxygen_aqueous(),
        use_mass_start=false) "Free dissolved O2 in water at 0degC" annotation (Placement(transformation(extent={{116,-34},{136,-14}})));
      Obsolete.Components.Substance water_0(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{126,-70},{146,-50}})));
      Obsolete.Sources.ExternalIdealGasSubstance CO2_g_0(
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
        TotalPressure=system.p_ambient,
        PartialPressure(displayUnit="mmHg") = 5332.8954966,
        Temperature=273.15) annotation (Placement(transformation(extent={{74,68},{94,88}})));
      Obsolete.Sources.ExternalIdealGasSubstance O2_g_0(
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        TotalPressure=system.p_ambient,
        PartialPressure(displayUnit="mmHg") = 12665.626804425,
        Temperature=273.15) annotation (Placement(transformation(extent={{112,68},{132,88}})));
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
</html>", revisions="<html>
<p><i>2015-2019</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        Diagram(coordinateSystem(extent={{-160,-100},{160,100}})));
    end GasSolubility_NIST;

    model GasSolubility "Dissolution of gases in liquids"
       extends Modelica.Icons.Example;

      Chemical.Obsolete.Components.Solution blood_plasma annotation (Placement(transformation(extent={{-100,-76},{-8,14}})));
                                          //(amountOfSolution_start=52.3)
      Chemical.Obsolete.Components.Solution red_cells annotation (Placement(transformation(extent={{8,-78},{102,14}})));
                                       //(amountOfSolution_start=39.7)
      Chemical.Obsolete.Components.GasSolubility CO2_dissolutionP annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
      Chemical.Obsolete.Components.Substance CO2_unbound_plasma(
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
        use_mass_start=false) "Free dissolved CO2 in blood plasma" annotation (Placement(transformation(extent={{-90,-24},{-70,-4}})));
      Chemical.Obsolete.Components.GasSolubility O2_dissolutionP annotation (Placement(transformation(extent={{-34,44},{-14,64}})));

      Chemical.Obsolete.Sources.ExternalIdealGasSubstance O2_g_n1(
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        PartialPressure=12665.626804425,
        TotalPressure=system.p_ambient) annotation (Placement(transformation(extent={{22,76},{42,96}})));
      Chemical.Obsolete.Components.Substance O2_unbound_plasma(
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        substanceData=Chemical.Obsolete.Substances.Oxygen_aqueous(),
        use_mass_start=false) "Free dissolved O2 in blood plasma" annotation (Placement(transformation(extent={{-50,-26},{-30,-6}})));
      Chemical.Obsolete.Components.GasSolubility CO2_dissolutionE annotation (Placement(transformation(extent={{36,44},{56,64}})));

      Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO2_g_n2(
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
        TotalPressure=system.p_ambient,
        PartialPressure(displayUnit="mmHg") = 5332.8954966) annotation (Placement(transformation(extent={{-58,78},{-38,98}})));

      Chemical.Obsolete.Components.Substance CO2_unbound_erythrocyte(
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
        use_mass_start=false) "Free dissolved CO2 in red cells" annotation (Placement(transformation(extent={{18,-32},{38,-12}})));

      Chemical.Obsolete.Components.GasSolubility O2_dissolutionE_NIST annotation (Placement(transformation(extent={{78,44},{98,64}})));
      Chemical.Obsolete.Components.Substance O2_unbound_erythrocyte_NIST(
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        substanceData=Chemical.Obsolete.Substances.Oxygen_aqueous(),
        use_mass_start=false) "Free dissolved O2 in red cells" annotation (Placement(transformation(extent={{58,-32},{78,-12}})));
      Obsolete.Components.Substance water_plasma(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=0.82)
        annotation (Placement(transformation(extent={{-40,-66},{-20,-46}})));
      Obsolete.Components.Substance water(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{68,-68},{88,-48}})));
      inner Modelica.Fluid.System system(p_ambient(displayUnit="mmHg")=
          101325.0144354, T_ambient=310.15)
        annotation (Placement(transformation(extent={{-10,-96},{10,-76}})));
      Obsolete.Components.Substance water_plasma1(mass_start=0.18, substanceData=Chemical.Obsolete.Interfaces.Incompressible.SubstanceData(MolarWeight=1/0.627))
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
</html>", revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end GasSolubility;

    model GasSolubility_blood "Dissolution of gases in liquids"
       extends Modelica.Icons.Example;

      Chemical.Obsolete.Components.Solution blood_plasma annotation (Placement(transformation(extent={{-100,-76},{-8,14}})));
                                          //(amountOfSolution_start=52.3)
      //(amountOfSolution_start=39.7)
      Chemical.Obsolete.Components.GasSolubility CO2_dissolutionP annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
      Chemical.Obsolete.Components.Substance CO2_unbound_plasma(
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
        use_mass_start=false) "Free dissolved CO2 in blood plasma" annotation (Placement(transformation(extent={{-90,-24},{-70,-4}})));
      Chemical.Obsolete.Components.GasSolubility O2_dissolutionP annotation (Placement(transformation(extent={{-34,44},{-14,64}})));

      Chemical.Obsolete.Sources.ExternalIdealGasSubstance O2_g_n1(
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        PartialPressure=12665.626804425,
        TotalPressure=system.p_ambient) annotation (Placement(transformation(extent={{22,76},{42,96}})));
      Chemical.Obsolete.Components.Substance O2_unbound_plasma(
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        substanceData=Chemical.Obsolete.Substances.Oxygen_aqueous(),
        use_mass_start=false) "Free dissolved O2 in blood plasma" annotation (Placement(transformation(extent={{-50,-26},{-30,-6}})));

      Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO2_g_n2(
        substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
        TotalPressure=system.p_ambient,
        PartialPressure(displayUnit="mmHg") = 5332.8954966) annotation (Placement(transformation(extent={{-58,78},{-38,98}})));

      Obsolete.Components.Substance water(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=0.82)
        annotation (Placement(transformation(extent={{-40,-66},{-20,-46}})));
      inner Modelica.Fluid.System system(p_ambient(displayUnit="mmHg")=
          101325.0144354, T_ambient=310.15)
        annotation (Placement(transformation(extent={{-10,-96},{10,-76}})));
      Obsolete.Components.Substance others(mass_start=0.18, substanceData=Chemical.Obsolete.Interfaces.Incompressible.SubstanceData(MolarWeight=1/0.627))
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
</html>", revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end GasSolubility_blood;

    model EnzymeKinetics "Basic enzyme kinetics"
      extends Modelica.Icons.Example;

      Chemical.Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
      Chemical.Obsolete.Components.Substance P(
        substanceData(DfG=-Modelica.Constants.R*298.15*50),
        use_mass_start=false,
        amountOfSubstance_start=1e-8) annotation (Placement(transformation(extent={{92,-12},{72,8}})));
      Chemical.Obsolete.Components.Substance S(use_mass_start=false, amountOfSubstance_start=1)
        annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));

      parameter Modelica.Units.SI.AmountOfSubstance tE=1e-6
        "Total amount of enzyme";
         parameter Real k_cat(unit="1/s", displayUnit="1/min")= 1
        "Forward rate of second reaction";
      constant Modelica.Units.SI.Concentration Km=0.1
        "Michaelis constant = substrate concentration at rate of half Vmax";

      parameter Modelica.Units.SI.MolarFlowRate Vmax=1e-5*k_cat
        "Maximal molar flow";

      Chemical.Obsolete.Components.Substance ES(
        substanceData(DfG=-Modelica.Constants.R*298.15*log(2/Km)),
        use_mass_start=false,
        amountOfSubstance_start=tE/2) annotation (Placement(transformation(extent={{-12,-10},{8,10}})));
      Chemical.Obsolete.Components.Substance E(use_mass_start=false, amountOfSubstance_start=tE/2)
        annotation (Placement(transformation(extent={{-10,38},{10,58}})));
      Chemical.Obsolete.Components.Reaction chemicalReaction(
        nS=2,
        KC=Vmax/(2*Modelica.Constants.R*298.15*log(2)),
        nP=1) annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));

      Chemical.Obsolete.Components.Reaction chemicalReaction1(
        nP=2,
        KC=Vmax/(2*Modelica.Constants.R*298.15*(50 - log(2))),
        nS=1) annotation (Placement(transformation(extent={{24,-10},{44,10}})));

      Obsolete.Components.Substance liquidWater(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
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
</html>",     info="<html>
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

      Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
      Chemical.Obsolete.Components.Substance Ag(
        substanceData=Chemical.Obsolete.Substances.Silver_solid(),
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{-72,-30},{-52,-10}})));
      Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{62,-50},{96,50}})));

      Chemical.Obsolete.Components.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-30,-60},{38,6}})));

      Chemical.Obsolete.Components.Substance Cl(
        substanceData=Chemical.Obsolete.Substances.Chloride_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=12.39) annotation (Placement(transformation(extent={{0,-26},{-20,-6}})));
      Chemical.Obsolete.Components.Substance AgCl(
        substanceData=Chemical.Obsolete.Substances.SilverChloride_solid(),
        use_mass_start=false,
        amountOfSubstance_start=1e-8) annotation (Placement(transformation(extent={{-76,4},{-56,24}})));
      Chemical.Obsolete.Sources.ExternalIdealGasSubstance H2(
        substanceData=Chemical.Obsolete.Substances.Hydrogen_gas(),
        PartialPressure=100000,
        TotalPressure=100000) annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Chemical.Obsolete.Components.Substance H(
        substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=12.39) annotation (Placement(transformation(extent={{8,-26},{28,-6}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Chemical.Obsolete.Components.Reaction electrodeReaction(
        nP=2,
        p={2,2},
        nS=1) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={52,6})));
      Chemical.Obsolete.Components.Reaction electrodeReaction1(nS=2, nP=2)
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-40,0})));

      Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                                 //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                  //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{84,-84},{104,-64}})));
      Obsolete.Components.Substance liquidWater(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
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

      Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{24,-76},{58,32}})));

      Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-80,-78},{-46,30}})));

      Chemical.Obsolete.Components.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-26,-80},{2,20}})));

      Chemical.Obsolete.Components.Substance Pb(
        substanceData=Chemical.Obsolete.Substances.Lead_solid(),
        use_mass_start=false,
        amountOfSubstance_start=50) annotation (Placement(transformation(extent={{50,-66},{30,-46}})));
      Chemical.Obsolete.Components.Substance HSO4(
        substanceData=Chemical.Obsolete.Substances.HydrogenSulfate_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{-2,-70},{-22,-50}})));
      Chemical.Obsolete.Components.Substance PbSO4_(
        amountOfSubstance_start(displayUnit="mol") = 0.001,
        substanceData=Chemical.Obsolete.Substances.LeadSulfate_solid(),
        use_mass_start=false) annotation (Placement(transformation(extent={{52,-30},{32,-10}})));
      Chemical.Obsolete.Components.Substance H(
        substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{-2,-42},{-22,-22}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-32,72},{-12,92}})));
      Chemical.Obsolete.Components.Reaction electrodeReaction(
        nP=2,
        nS=4,
        s={1,1,3,2},
        p={1,2}) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-36,-14})));
      Chemical.Obsolete.Components.Reaction electrodeReaction1(
        nS=2,
        nP=3,
        p={1,1,2}) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={14,-16})));

      Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{50,2},{30,22}})));
      Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{-72,-38},{-52,-18}})));
      Chemical.Obsolete.Components.Substance PbO2(
        substanceData=Chemical.Obsolete.Substances.LeadDioxide_solid(),
        use_mass_start=false,
        amountOfSubstance_start=50) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-60,-58})));
      Obsolete.Components.Substance H2O(mass_start(displayUnit="g") = 0.114, substanceData=Chemical.Obsolete.Substances.Water_liquid())
        annotation (Placement(transformation(extent={{-2,-8},{-22,12}})));
      Chemical.Obsolete.Components.Substance PbSO4(
        amountOfSubstance_start=0.001,
        substanceData=Chemical.Obsolete.Substances.LeadSulfate_solid(),
        use_mass_start=false) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-60,6})));

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
</html>",   info="<html>
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

        Chemical.Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-72,2},{76,96}})));
        Chemical.Obsolete.Components.Solution solution1 annotation (Placement(transformation(extent={{-76,-98},{72,-4}})));
        Chemical.Obsolete.Components.Substance H3O(
          substanceData=Chemical.Obsolete.Substances.Hydronium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={30,70})));
        Chemical.Obsolete.Components.Substance OH(
          substanceData=Chemical.Obsolete.Substances.Hydroxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={30,26})));
        Chemical.Obsolete.Components.Substance H2O(mass_start=1, substanceData=Chemical.Obsolete.Substances.Water_liquid())
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-30,46})));
        Chemical.Obsolete.Components.Reaction waterDissociation(
          KC=1e-3,
          nS=1,
          nP=2,
          s={2}) annotation (Placement(transformation(extent={{-12,36},{8,56}})));
              Real pH, pH3O;
        Chemical.Obsolete.Components.Substance H_(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={28,-30})));
        Chemical.Obsolete.Components.Substance OH_(
          substanceData=Chemical.Obsolete.Substances.Hydroxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={28,-76})));
        Chemical.Obsolete.Components.Substance H2O_(mass_start=1, substanceData=Chemical.Obsolete.Substances.Water_liquid())
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-32,-56})));
        Chemical.Obsolete.Components.Reaction waterDissociation_(
          KC=1e-3,
          nS=1,
          nP=2) annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));

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
</html>",      revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=1));
      end WaterSelfIonization;

      model CarbonDioxideInWater "CO2 as alone acid-base buffer"
        import Chemical;
          extends Modelica.Icons.Example;
          parameter Real KC = 1e-3;

        Chemical.Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,46}})));
        Chemical.Obsolete.Components.Substance HCO3(
          amountOfSubstance_start(displayUnit="mmol") = 1e-08,
          substanceData=Chemical.Obsolete.Substances.Bicarbonate_aqueous(),
          use_mass_start=false) annotation (Placement(transformation(extent={{-16,-4},{4,16}})));
        Chemical.Obsolete.Components.Reaction HendersonHasselbalch(
          KC=1e-4*KC,
          nP=2,
          nS=2,
          useKineticsInput=false) "K=10^(-6.103 + 3), dH=7.3 kJ/mol" annotation (Placement(transformation(extent={{-48,-6},{-28,14}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO2_gas(
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
          PartialPressure(displayUnit="mmHg") = 5332.8954966,
          TotalPressure(displayUnit="mmHg") = 101325.0144354)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-60,86})));
        Chemical.Obsolete.Components.Substance H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={10,-30})));
        Chemical.Obsolete.Components.GasSolubility gasSolubility(KC=KC) annotation (Placement(transformation(extent={{-70,36},{-50,56}})));
                                              /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
        Chemical.Obsolete.Components.Substance CO2_liquid(
          amountOfSubstance_start(displayUnit="mmol") = 0.001,
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
          use_mass_start=false) annotation (Placement(transformation(extent={{-82,-6},{-62,14}})));
        Real pH;

        Chemical.Obsolete.Components.Substance liquidWater(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-76,-50},{-56,-30}})));
        inner Modelica.Fluid.System system(T_ambient=310.15)
          annotation (Placement(transformation(extent={{48,64},{68,84}})));
        Chemical.Obsolete.Components.Substance OH(
          substanceData=Chemical.Obsolete.Substances.Hydroxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={12,-72})));
        Chemical.Obsolete.Components.Reaction waterDissociation(
          KC=KC,
          nP=2,
          nS=1,
          useKineticsInput=false) "water dissociation" annotation (Placement(transformation(extent={{-44,-68},{-24,-48}})));
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
</html>",      revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.02));
      end CarbonDioxideInWater;

      model Phosphate
        import Chemical;
          extends Modelica.Icons.Example;
          parameter Real KC = 1e-3;
        Chemical.Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-98,-100},{100,100}})));

        Chemical.Obsolete.Components.Substance H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=10^(-7.4)) "hydrogen ions activity" annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={28,-14})));

        Chemical.Obsolete.Components.Substance H3PO4(
          substanceData=Chemical.Obsolete.Substances.PhosphoricAcid_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-08) annotation (Placement(transformation(extent={{-90,-58},{-70,-38}})));
        Chemical.Obsolete.Components.Substance H2PO4(
          substanceData=Chemical.Obsolete.Substances.DihydrogenPhosphate_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.0005) annotation (Placement(transformation(extent={{-40,-58},{-20,-38}})));
        Chemical.Obsolete.Components.Substance HPO4(
          substanceData=Chemical.Obsolete.Substances.HydrogenPhosphate_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.0006) annotation (Placement(transformation(extent={{16,-58},{36,-38}})));
        Chemical.Obsolete.Components.Substance PO4(
          substanceData=Chemical.Obsolete.Substances.Phosphate_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-08) annotation (Placement(transformation(extent={{92,-58},{72,-38}})));

        Chemical.Obsolete.Components.Reaction chemicalReaction(
          KC=KC,
          nS=1,
          nP=2) "10^(-1.915 + 3)" annotation (Placement(transformation(extent={{-66,-58},{-46,-38}})));
        Chemical.Obsolete.Components.Reaction chemicalReaction1(
          KC=KC,
          nS=1,
          nP=2) "10^(-6.66 + 3)" annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
        Chemical.Obsolete.Components.Reaction chemicalReaction2(
          KC=KC,
          nS=1,
          nP=2) "10^(-11.78 + 3)" annotation (Placement(transformation(extent={{44,-58},{64,-38}})));

        Chemical.Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={86,14})));
        Real pH "acidity";
        Chemical.Obsolete.Components.Substance OH(
          substanceData=Chemical.Obsolete.Substances.Hydroxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=10^(-14 + 7.4)) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={28,32})));
        Chemical.Obsolete.Components.Reaction reaction(
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
</html>",      revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.05));
      end Phosphate;

      model CarbonDioxideInBlood
        import Chemical;
          extends Modelica.Icons.Example;

        parameter Real KC=1e-3;
        //e-6 "Slow down factor";
        Chemical.Obsolete.Components.Solution blood_erythrocytes(ElectricGround=false, temperature_start=310.15)
          annotation (Placement(transformation(extent={{-100,-98},{100,-38}})));
        Chemical.Obsolete.Components.Solution blood_plasma(temperature_start=310.15) annotation (Placement(transformation(extent={{-100,4},{100,56}})));

        Chemical.Obsolete.Components.Substance HCO3(
          substanceData=Chemical.Obsolete.Substances.Bicarbonate_blood(),
          use_mass_start=false,
          amountOfSubstance_start=0.024) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={18,24})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO2_gas(
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
          TotalPressure(displayUnit="mmHg") = 101325.0144354,
          PartialPressure(displayUnit="mmHg") = 5332.8954966,
          usePartialPressureInput=true,
          Temperature=310.15) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,84})));
        Chemical.Obsolete.Components.GasSolubility gasSolubility(KC=KC) annotation (Placement(transformation(extent={{-94,48},{-74,68}})));

        Chemical.Obsolete.Components.Substance CO2(
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00148) "Free dissolved CO2 in plasma" annotation (Placement(transformation(extent={{-88,28},{-68,48}})));
        Chemical.Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=51.6159/55.508)
          annotation (Placement(transformation(extent={{-60,12},{-40,32}})));
        Chemical.Obsolete.Components.Substance HCO3_E(
          substanceData=Chemical.Obsolete.Substances.Bicarbonate_blood(),
          use_mass_start=false,
          amountOfSubstance_start=0.0116) annotation (Placement(transformation(extent={{28,-60},{8,-40}})));
        Chemical.Obsolete.Components.Reaction HendersonHasselbalch1(
          nP=2,
          nS=2,
          KC=KC) "K=10^(-6.103 + 3), dH=7.3 kJ/mol" annotation (Placement(transformation(extent={{-26,-68},{-6,-48}})));
        Chemical.Obsolete.Components.Substance CO2_E(
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.0011) "Free dissolved CO2 in erythrocyte" annotation (Placement(transformation(extent={{-90,-82},{-70,-62}})));
        Chemical.Obsolete.Components.Substance H2O_E(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=38.4008/55.508)
          annotation (Placement(transformation(extent={{-60,-62},{-40,-42}})));
        Chemical.Obsolete.Components.Substance Cl_E(
          substanceData=Chemical.Obsolete.Substances.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.0499) annotation (Placement(transformation(extent={{68,-60},{48,-40}})));
        Chemical.Obsolete.Components.Substance Cl(
          substanceData=Chemical.Obsolete.Substances.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.103) annotation (Placement(transformation(extent={{68,20},{48,40}})));

        Real pH_e, pH_p;

        Chemical.Obsolete.Components.Membrane aquaporin(KC=KC)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-34,-16})));
        Chemical.Obsolete.Components.Membrane Band3_HCO3(KC=KC)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={4,-16})));
        Chemical.Obsolete.Components.Membrane Band3_Cl(useKineticsInput=false, KC=KC)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={46,-16})));
        Chemical.Obsolete.Sources.Buffer H_E(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          BufferValue=0.063,
          a_start=10^(-7.2)) annotation (Placement(transformation(extent={{48,-84},{30,-66}})));
        Modelica.Blocks.Sources.ContinuousClock clock(offset=5000)
          annotation (Placement(transformation(extent={{-54,62},{-34,82}})));
        Chemical.Obsolete.Components.Substance others_E(
          substanceData=Chemical.Obsolete.Interfaces.Incompressible.SubstanceData(
                    density=(1.045 - 0.695523)*1000/(1 - 0.697583),
                    References={"erythrocyte intracellular fluid density 1045kg/m3"},
                    MolarWeight=(1.045 - 0.695523)/(38.7*(1 - 0.994648) - 0.0499 - 0.0116 - 0.00123)),
          use_mass_start=false,
          amountOfSubstance_start=0.1444) annotation (Placement(transformation(extent={{68,-88},{88,-68}})));
        Chemical.Obsolete.Components.Substance others_P(
          substanceData=Chemical.Obsolete.Interfaces.Incompressible.SubstanceData(
                    References={"to reach plasma density 1024 kg/m3 and plasma volume 1 liter"},
                    density=(1.024 - 0.933373)*1000/(1 - 0.936137),
                    MolarWeight=(1.024 - 0.933373)/(51.8*(1 - 0.994648) - 0.103 - 0.024 - 0.0017)),
          use_mass_start=false,
          amountOfSubstance_start=0.1487) annotation (Placement(transformation(extent={{70,14},{90,34}})));
        Chemical.Obsolete.Components.Diffusion diffusion
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-66,-16})));
        Chemical.Obsolete.Sources.Buffer H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          BufferValue=0.0077,
          a_start=10^(-7.4)) "buffer value 7.7 mmol/L for plasma is from (O. Siggaard-Andersen 1995)"
          annotation (Placement(transformation(extent={{38,38},{20,56}})));
        Chemical.Obsolete.Components.Reaction HendersonHasselbalch2(
          nP=2,
          nS=2,
          KC=(1e-10)*KC) "K=10^(-6.103 + 3), dH=7.3 kJ/mol" annotation (Placement(transformation(extent={{-26,26},{-6,46}})));
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
</html>",      revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(
          StopTime=3),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));
      end CarbonDioxideInBlood;

      model AcidBaseBufferTest
          extends Modelica.Icons.Example;

        Chemical.Obsolete.Sources.Buffer buffer(
          substanceData(z=1.045),
          a_start=10^(-7.2),
          BufferValue=3) annotation (Placement(transformation(extent={{-50,4},{-30,24}})));
        Chemical.Obsolete.Components.Solution simpleSolution annotation (Placement(transformation(extent={{-104,-100},{96,100}})));
        Chemical.Obsolete.Sources.ExternalMoleFraction externalMoleFraction(substanceData=Chemical.Obsolete.Substances.Proton_aqueous(), MoleFraction=10^(-7.1))
          annotation (Placement(transformation(extent={{0,-46},{20,-26}})));
        Obsolete.Components.Substance liquidWater(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
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
          Chemical.Obsolete.Components.Solution blood_erythrocytes(ElectricGround=false) annotation (Placement(transformation(extent={{-180,-100},{180,-10}})));
          Chemical.Obsolete.Components.Solution blood_plasma annotation (Placement(transformation(extent={{-180,12},{180,100}})));

          Chemical.Obsolete.Components.Substance HCO3(
            substanceData=Chemical.Obsolete.Substances.Bicarbonate_blood(),
            use_mass_start=false,
            amountOfSubstance_start=0.024) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-18,30})));

          Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=51.8*0.994648/55.508)
            annotation (Placement(transformation(extent={{-146,44},{-166,64}})));
          Chemical.Obsolete.Components.Substance HCO3_E(
            substanceData=Chemical.Obsolete.Substances.Bicarbonate_blood(),
            use_mass_start=false,
            amountOfSubstance_start=0.0116) annotation (Placement(transformation(extent={{-28,-38},{-8,-18}})));
          Obsolete.Components.Substance H2O_E(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=38.7*0.994648/55.508)
            annotation (Placement(transformation(extent={{-144,-38},{-164,-18}})));
          Chemical.Obsolete.Components.Substance Cl_E(
            substanceData=Chemical.Obsolete.Substances.Chloride_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.0499) annotation (Placement(transformation(extent={{-4,-38},{16,-18}})));
          Chemical.Obsolete.Components.Substance Cl(
            substanceData=Chemical.Obsolete.Substances.Chloride_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.103) annotation (Placement(transformation(extent={{-4,20},{16,40}})));
          Chemical.Obsolete.Components.Substance albumin(
            substanceData(
              MolarWeight=66.463,
              z=-17,
              density=1080),
            use_mass_start=false,
            amountOfSubstance_start=0.0007) annotation (Placement(transformation(extent={{112,76},{92,96}})));
          Real pH_e,pH_p;

          Chemical.Obsolete.Components.Membrane Aquapirin(KC=KC)
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-168,0})));
          Chemical.Obsolete.Components.Membrane Band3(KC=KC)
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-6,0})));
          Chemical.Obsolete.Components.Membrane Band3_(useKineticsInput=false, KC=KC)
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={18,0})));
          Chemical.Obsolete.Components.Substance permeableUncharged(use_mass_start=false, amountOfSubstance_start=0.0118)
            annotation (Placement(transformation(extent={{166,20},{146,40}})));
          Chemical.Obsolete.Components.Substance permeableUncharged_E(
            substanceData(MolarWeight=0.1),
            use_mass_start=false,
            amountOfSubstance_start=0.00903) annotation (Placement(transformation(extent={{164,-38},{144,-18}})));
          Chemical.Obsolete.Components.Substance chargedImpermeable_E(
            substanceData(MolarWeight=1),
            use_mass_start=false,
            amountOfSubstance_start=0.0165) annotation (Placement(transformation(extent={{144,-62},{164,-42}})));
          Chemical.Obsolete.Components.Membrane leak(useKineticsInput=false, KC=KC)
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={140,0})));
          Chemical.Obsolete.Components.Substance Lac_E(
            substanceData=Chemical.Obsolete.Substances.Chloride_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.00062) annotation (Placement(transformation(extent={{56,-38},{76,-18}})));
          Chemical.Obsolete.Components.Substance Lac(
            substanceData=Chemical.Obsolete.Substances.Chloride_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.00131) annotation (Placement(transformation(extent={{56,20},{76,40}})));
          Chemical.Obsolete.Components.Membrane MCT_(useKineticsInput=false, KC=KC) "Monocarboxylate transporters"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={78,0})));
          Chemical.Obsolete.Components.Substance H_E(
            substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=10^(-7.2)) "H+" annotation (Placement(transformation(extent={{30,-38},{50,-18}})));
          Chemical.Obsolete.Components.Substance H(
            substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=10^(-7.4)) "H+ in plasma" annotation (Placement(transformation(extent={{30,20},{50,40}})));
          Chemical.Obsolete.Components.Membrane MCT(useKineticsInput=false, KC=KC) "Monocarboxylate transporters"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={52,0})));
          Chemical.Obsolete.Components.Substance CO2(
            substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.00167) "free dissolved unbound CO2" annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
          Chemical.Obsolete.Components.Substance CO2_E(
            substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.00125) "free dissolved unbound CO2" annotation (Placement(transformation(extent={{-58,-38},{-38,-18}})));
          Chemical.Obsolete.Components.Membrane freeCO2(KC=KC)
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-38,2})));
          Chemical.Obsolete.Components.Substance O2(
            substanceData=Chemical.Obsolete.Substances.Oxygen_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.000167) "free dissolved undound oxygen" annotation (Placement(transformation(extent={{96,20},{116,40}})));
          Chemical.Obsolete.Components.Membrane freeO2(KC=KC)
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={118,0})));
          Chemical.Obsolete.Components.Substance O2_E(
            substanceData=Chemical.Obsolete.Substances.Oxygen_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.000125) "free dissolved undound O2" annotation (Placement(transformation(extent={{96,-38},{116,-18}})));
          Chemical.Obsolete.Components.Substance K(
            substanceData=Chemical.Obsolete.Substances.Potassium_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.004) annotation (Placement(transformation(extent={{-100,20},{-120,40}})));
          Chemical.Obsolete.Components.Substance Na(
            substanceData=Chemical.Obsolete.Substances.Sodium_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.138) annotation (Placement(transformation(extent={{-124,20},{-144,40}})));
          Chemical.Obsolete.Components.Substance Na_E(
            substanceData=Chemical.Obsolete.Substances.Sodium_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.007) annotation (Placement(transformation(extent={{-118,-38},{-138,-18}})));
          Chemical.Obsolete.Components.Substance K_E(
            substanceData=Chemical.Obsolete.Substances.Potassium_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.096) annotation (Placement(transformation(extent={{-112,-38},{-92,-18}})));
          Chemical.Obsolete.Components.Substance H2PO4_E(
            substanceData=Chemical.Obsolete.Substances.DihydrogenPhosphate_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.000175) annotation (Placement(transformation(extent={{-84,-38},{-64,-18}})));
          Chemical.Obsolete.Components.Substance ADP_E(
            substanceData(z=-3),
            use_mass_start=false,
            amountOfSubstance_start=9.6e-05) annotation (Placement(transformation(extent={{-114,-62},{-94,-42}})));
          Chemical.Obsolete.Components.Substance ATP_E(
            substanceData(
              z=-4,
              DfH=16700,
              DfG=30500,
              References={"http://www.wiley.com/college/pratt/0471393878/student/review/thermodynamics/7_relationship.html"}),
            use_mass_start=false,
            amountOfSubstance_start=0.00128) annotation (Placement(transformation(extent={{-146,-62},{-166,-42}})));

          Chemical.Obsolete.Components.Substance HPO4_E(
            substanceData=Chemical.Obsolete.Substances.HydrogenPhosphate_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.000495) annotation (Placement(transformation(extent={{-84,-62},{-64,-42}})));
          Chemical.Obsolete.Components.Substance globulins(
            substanceData(
              MolarWeight=34,
              z=-2.43,
              density=1080),
            use_mass_start=false,
            amountOfSubstance_start=0.00082) annotation (Placement(transformation(extent={{150,76},{130,96}})));
          Chemical.Obsolete.Components.Substance Ca(
            substanceData=Chemical.Obsolete.Substances.Calcium_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.00175) "Ca2+" annotation (Placement(transformation(extent={{-78,20},{-98,40}})));
          Chemical.Obsolete.Components.Substance Mg(
            substanceData=Chemical.Obsolete.Substances.Magnesium_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.00108) "Mg2+" annotation (Placement(transformation(extent={{-112,-84},{-92,-64}})));
          Chemical.Obsolete.Components.Substance DPG(
            substanceData(
              MolarWeight=0.266,
              z=-2.2,
              density=1000),
            use_mass_start=false,
            amountOfSubstance_start=0.0051) annotation (Placement(transformation(extent={{128,-94},{108,-74}})));
          Chemical.Obsolete.Components.Substance GSH(
            substanceData(
              MolarWeight=0.2,
              z=-1,
              density=1000),
            use_mass_start=false,
            amountOfSubstance_start=0.00223) annotation (Placement(transformation(extent={{164,-94},{144,-74}})));
          Obsolete.Components.Reaction HendersonHasselbalch(
            nP=2,
            nS=2,
            useKineticsInput=false) "K=10^(-6.103 + 3), dH=7.3 kJ/mol" annotation (Placement(transformation(extent={{-34,64},{-14,44}})));
          Obsolete.Sources.Buffer Hemoglobin(
            substanceData(z=1.045),
            a_start=10^(-7.2),
            BufferValue=3) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={82,-76})));
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
</html>",        revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),  experiment(StopTime=1e-008),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},
                    {180,100}}), graphics));
        end RedCellMembrane;

        model NaKATPase
          Obsolete.Components.Solution ICF(ElectricGround=false) annotation (Placement(transformation(extent={{-100,-100},{-20,98}})));
          Obsolete.Components.Solution ECF annotation (Placement(transformation(extent={{28,-100},{100,100}})));
          Obsolete.Components.Substance water_ICF(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
            annotation (Placement(transformation(extent={{-82,70},{-62,90}})));
          Obsolete.Components.Substance water_EFC(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
            annotation (Placement(transformation(extent={{58,56},{38,76}})));
          Obsolete.Components.Reaction NaKATPase(
            nS=4,
            nP=4,
            s={3,2,1,1},
            p={3,2,1,1}) annotation (Placement(transformation(
                extent={{-10,10},{10,-10}},
                rotation=270,
                origin={0,8})));
          Obsolete.Components.Substance Na_ICF(
            amountOfSubstance_start(displayUnit="mmol") = 0.01,
            substanceData=Chemical.Obsolete.Substances.Sodium_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-82,38},{-62,58}})));
          Obsolete.Components.Substance K_ICF(
            amountOfSubstance_start(displayUnit="mmol") = 0.14,
            substanceData=Chemical.Obsolete.Substances.Potassium_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-82,10},{-62,30}})));
          Obsolete.Components.Substance ATP4_ICF(
            amountOfSubstance_start(displayUnit="mmol") = 0.0038,
            substanceData=Chemical.Obsolete.Substances.ATP4_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-82,-16},{-62,4}})));
          Obsolete.Components.Substance Na_ECF(
            amountOfSubstance_start(displayUnit="mmol") = 0.14,
            substanceData=Chemical.Obsolete.Substances.Sodium_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{58,-26},{38,-6}})));
          Obsolete.Components.Substance K_ECF(
            amountOfSubstance_start(displayUnit="mmol") = 0.005,
            substanceData=Chemical.Obsolete.Substances.Potassium_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{58,28},{38,48}})));
          Obsolete.Components.Substance ADP3(
            substanceData=Chemical.Obsolete.Substances.ADP3_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.005/150) annotation (Placement(transformation(extent={{-82,-42},{-62,-22}})));
          Obsolete.Components.Substance H2PO4_ICF(
            amountOfSubstance_start(displayUnit="mmol") = 0.036,
            substanceData=Chemical.Obsolete.Substances.DihydrogenPhosphate_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-82,-68},{-62,-48}})));
          inner Modelica.Fluid.System system(T_ambient=310.15)
            annotation (Placement(transformation(extent={{-6,-90},{14,-70}})));
          Obsolete.Components.Substance Cl_ICF(
            amountOfSubstance_start(displayUnit="mmol") = 0.0987,
            substanceData=Chemical.Obsolete.Substances.Chloride_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-82,-92},{-62,-72}})));
          Obsolete.Components.Substance Cl_ECF(
            amountOfSubstance_start(displayUnit="mmol") = 0.145,
            substanceData=Chemical.Obsolete.Substances.Chloride_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{50,-72},{70,-52}})));
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

        Chemical.Obsolete.Components.Solution solution(redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible)
          annotation (Placement(transformation(extent={{-104,-100},{96,100}})));

        constant Integer n=218 "Number of weak acid group in albumin molecule";
        constant Real pKAs[n]=cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
          "acid dissociation constants";
        constant Real K[n]=fill(10.0, n) .^ (-pKAs);
        constant Real DfG[n]= Modelica.Constants.R*(298.15)*log(K);

        Chemical.Obsolete.Components.Substance A[n](
          substanceData(each z=-1),
          each use_mass_start=false,
          each amountOfSubstance_start=0.00033) "deprotonated acid groups" annotation (Placement(transformation(extent={{26,-16},{6,4}})));
        Chemical.Obsolete.Components.Reaction react[n](
          each KC=1e-9,
          each nS=1,
          each nP=2) annotation (Placement(transformation(extent={{-44,-2},{-24,18}})));

        Chemical.Obsolete.Components.Substance HA[n](
          substanceData(DfG=DfG),
          each use_mass_start=false,
          each amountOfSubstance_start=0.00033) "protonated acid groups" annotation (Placement(transformation(extent={{-78,-2},{-58,18}})));

        Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,-68})));
        Obsolete.Sources.ExternalMoleFraction H(substanceData=Chemical.Obsolete.Substances.Proton_aqueous(), MoleFraction=10^(-7.4))
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
</html>",       info="<html>
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

        Chemical.Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-72,-102},{94,124}})));

        Chemical.Obsolete.Components.Substance oxygen_unbound(
          substanceData(DfG=GO2aq),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-62,-46},{-42,-26}})));

        Chemical.Obsolete.Components.Substance T0(
          substanceData(DfG=GT0),
          use_mass_start=false,
          amountOfSubstance_start=(THb)) annotation (Placement(transformation(extent={{34,78},{54,98}})));

        Chemical.Obsolete.Components.Substance T1(
          substanceData(DfG=GT1),
          use_mass_start=false,
          amountOfSubstance_start=(THb*1e-4)) annotation (Placement(transformation(extent={{34,36},{54,56}})));

        Chemical.Obsolete.Components.Substance T2(
          substanceData(DfG=GT2),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-8) annotation (Placement(transformation(extent={{34,-10},{54,10}})));

        Chemical.Obsolete.Components.Substance R1(
          substanceData(DfG=GR1),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-8) annotation (Placement(transformation(extent={{-20,36},{0,56}})));

        Chemical.Obsolete.Components.Substance R2(
          substanceData(DfG=GR2),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-10) annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

        Chemical.Obsolete.Components.Substance T3(
          substanceData(DfG=GT3),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-12) annotation (Placement(transformation(extent={{34,-54},{54,-34}})));

        Chemical.Obsolete.Components.Substance R3(
          substanceData(DfG=GR3),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-12) annotation (Placement(transformation(extent={{-20,-54},{0,-34}})));

        Chemical.Obsolete.Components.Substance T4(
          substanceData(DfG=GT4),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-17) annotation (Placement(transformation(extent={{34,-92},{54,-72}})));

        Chemical.Obsolete.Components.Substance R4(
          substanceData(DfG=GR4),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-14) annotation (Placement(transformation(extent={{-20,-92},{0,-72}})));

        Chemical.Obsolete.Components.Substance R0(
          substanceData(DfG=GR0),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-7) annotation (Placement(transformation(extent={{-20,78},{0,98}})));

        Chemical.Obsolete.Components.Reaction quaternaryForm(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{4,78},{24,98}})));
        Chemical.Obsolete.Components.Reaction oxyR1(
          nS=1,
          nP=2,
          KC=KC) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-8,64})));
        Chemical.Obsolete.Components.Reaction oxyT1(
          nS=1,
          nP=2,
          KC=KC) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,64})));
        Chemical.Obsolete.Components.Reaction oxyR2(
          nS=1,
          nP=2,
          KC=KC) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,22})));
        Chemical.Obsolete.Components.Reaction oxyR3(
          nS=1,
          nP=2,
          KC=KC) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-24})));
        Chemical.Obsolete.Components.Reaction oxyR4(
          nS=1,
          nP=2,
          KC=KC) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-66})));
        Chemical.Obsolete.Components.Reaction oxyT2(
          nS=1,
          nP=2,
          KC=KC) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,22})));
        Chemical.Obsolete.Components.Reaction oxyT3(
          nS=1,
          nP=2,
          KC=KC) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-24})));
        Chemical.Obsolete.Components.Reaction oxyT4(
          nS=1,
          nP=2,
          KC=KC) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-66})));
        Chemical.Obsolete.Components.Reaction quaternaryForm1(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{8,36},{28,56}})));
        Chemical.Obsolete.Components.Reaction quaternaryForm2(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{8,-10},{28,10}})));
        Chemical.Obsolete.Components.Reaction quaternaryForm3(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{8,-54},{28,-34}})));
        Chemical.Obsolete.Components.Reaction quaternaryForm4(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{10,-92},{30,-72}})));

        Modelica.Blocks.Sources.ContinuousClock clock(offset=10)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,62})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance O2_in_air(
          usePartialPressureInput=true,
          TotalPressure(displayUnit="kPa") = 101325.0144354,
          substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
          PartialPressure(displayUnit="kPa") = 3733)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,22})));

        Chemical.Obsolete.Components.GasSolubility gasSolubility(KC=KC) annotation (Placement(transformation(extent={{-94,-16},{-74,4}})));

        Real sO2;
        Obsolete.Components.Substance substance(substanceData=Obsolete.Substances.Water_liquid(), mass_start=1)
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
</html>",   revisions="<html>
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
        Chemical.Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,42}})));

        Chemical.Obsolete.Components.Reaction quaternaryForm(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{12,-58},{32,-38}})));
        Chemical.Obsolete.Components.Speciation R0_in_R(NumberOfSubunits=4) annotation (Placement(transformation(extent={{-46,-48},{-26,-28}})));
         // AmountOfSubstance_start=4e-11)
        Chemical.Obsolete.Components.Speciation T0_in_T(NumberOfSubunits=4) annotation (Placement(transformation(extent={{76,-48},{56,-28}})));
         // AmountOfSubstance_start=totalAmountOfHemoglobin)
        Chemical.Obsolete.Components.Substance OxyRHm[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_O2 + RT*log(KRx) + DfG_tR/4),
          each use_mass_start=false,
          each amountOfSubstance_start=5.88e-9) "Oxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-96,-18},{-76,2}})));

        Chemical.Obsolete.Components.Reaction oxygenation_R[4](
          each nS=1,
          each nP=2,
          each KC=KC) annotation (Placement(transformation(extent={{-68,-18},{-48,2}})));
        Chemical.Obsolete.Components.Substance DeoxyRHm[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_tR/4),
          each use_mass_start=false,
          each amountOfSubstance_start=1.58e-7) "Deoxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-8,-20},{-28,0}})));

        Chemical.Obsolete.Components.Substance OxyTHm[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_O2 + RT*log(KTx) + DfG_tT/4),
          each use_mass_start=false,
          each amountOfSubstance_start=1e-4) "Oxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{14,-18},{34,2}})));

        Chemical.Obsolete.Components.Reaction oxygenation_T[4](
          each nS=1,
          each nP=2,
          each KC=KC) annotation (Placement(transformation(extent={{42,-18},{62,2}})));
        Chemical.Obsolete.Components.Substance DeoxyTHm[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_tT/4),
          each use_mass_start=false,
          each amountOfSubstance_start=THb - 1e-4 - 1.58e-7 - 5.88e-9) "Deoxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{96,-20},{76,0}})));

        Chemical.Obsolete.Components.Substance oxygen_unbound(
          substanceData(DfG=DfG_O2),
          use_mass_start=false,
          amountOfSubstance_start=2e-9) annotation (Placement(transformation(extent={{-2,6},{18,26}})));
        Modelica.Blocks.Sources.ContinuousClock clock(offset=1) annotation (
           Placement(transformation(extent={{-40,74},{-20,94}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance oxygen_in_air(usePartialPressureInput=true, substanceData=Chemical.Obsolete.Substances.Oxygen_gas())
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={8,68})));
        Chemical.Obsolete.Components.GasSolubility partialPressure1(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={8,40})));

        Real sO2 "Hemoglobin oxygen saturation";
        Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
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
</html>", info="<html>
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

        Chemical.Obsolete.Components.Speciation speciation(NumberOfSubunits=N) annotation (Placement(transformation(extent={{-18,-72},{2,-52}})));
         // AmountOfSubstance_start=4e-11)
        // AmountOfSubstance_start=totalAmountOfHemoglobin)
        Chemical.Obsolete.Components.Substance OxyHm[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_O2 + RT*log(Ko25) + DfG_selectedForm/N, DfH=DfH_O2 - Ho + DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=initialO2*initialHb/(Ko37 + initialO2)) "Oxygenated subunit"
          annotation (Placement(transformation(extent={{-88,14},{-68,34}})));

        Chemical.Obsolete.Components.Reaction o[4](
          each nS=1,
          each nP=2,
          each KC=KC) annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
        Chemical.Obsolete.Components.Substance DeoxyHm[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=Ko37*initialHb/(Ko37 + initialO2)) "Deoxygenated subunit"
          annotation (Placement(transformation(extent={{-8,12},{-28,32}})));

        Chemical.Obsolete.Interfaces.SolutionPort solution
          annotation (Placement(transformation(extent={{-50,-82},{-32,-62}}), iconTransformation(extent={{-50,-90},{-30,-70}})));
        Chemical.Obsolete.Interfaces.SubstancePort_b O2
          annotation (Placement(transformation(extent={{-28,32},{-8,52}}), iconTransformation(extent={{-90,70},{-70,90}})));
        Chemical.Obsolete.Interfaces.SubstancePort_a selectedForm
          annotation (Placement(transformation(extent={{26,-82},{46,-62}}), iconTransformation(extent={{30,-90},{50,-70}})));
        Chemical.Obsolete.Components.Substance HmAH[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=RT*log(Kh25) + DfG_selectedForm/N, DfH=-Hh + DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=initialH*initialHb/(Kh37 + initialH)) "Protonated h site of subunit in quaternary structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{10,12},{30,32}})));

        Chemical.Obsolete.Components.Reaction h[4](
          each nS=1,
          each nP=2,
          each KC=KC) annotation (Placement(transformation(extent={{36,32},{56,12}})));
        Chemical.Obsolete.Components.Substance HmA[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=Kh37*initialHb/(Kh37 + initialH)) "Deprotonated h site of subunit in quaternary structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{92,14},{72,34}})));

        Chemical.Obsolete.Components.Substance HmNH3[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=RT*log(Kz25) + DfG_selectedForm/N, DfH=-Hz + DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=(initialH^2)*initialHb/(initialH^2 + initialH*Kz37 + Kz37*Kc37*initialCO2))
          "Protonated z site of subunit in quaternary structure of hemoglobin tetramer" annotation (Placement(transformation(extent={{-84,-42},{-64,-22}})));

        Chemical.Obsolete.Components.Reaction z[4](
          each nS=1,
          each nP=2,
          each KC=KC) annotation (Placement(transformation(extent={{-54,-42},{-34,-22}})));
        Chemical.Obsolete.Components.Substance HmNH2[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=initialH*Kz37*initialHb/(initialH^2 + initialH*Kz37 + Kz37*Kc37*initialCO2))
          "Deprotonated z site of subunit in quaternary structure of hemoglobin tetramer" annotation (Placement(transformation(extent={{12,-44},{-8,-24}})));

        Chemical.Obsolete.Components.Reaction c[4](
          each nP=2,
          each KC=KC,
          each nS=2) annotation (Placement(transformation(extent={{20,-42},{40,-22}})));
        Chemical.Obsolete.Interfaces.SubstancePort_b CO2
          annotation (Placement(transformation(extent={{-8,-26},{12,-6}}), iconTransformation(extent={{10,70},{30,90}})));
        Chemical.Obsolete.Components.Substance HmNHCOO[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_CO2 - RT*log(Kc25) + DfG_selectedForm/N, DfH=DfH_CO2 + Hc + DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=(Kz37*Kc37*initialCO2)*initialHb/(initialH^2 + initialH*Kz37 + Kz37*Kc37*initialCO2))
          "Carboxylated c site of subunit in quaternary structure of hemoglobin tetramer" annotation (Placement(transformation(extent={{70,-44},{50,-24}})));

        Obsolete.Interfaces.SubstancePort_b H annotation (Placement(transformation(extent={{50,-12},{70,8}}), iconTransformation(extent={{70,70},{90,90}})));
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
</html>", info="<html>
<p>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
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

        Chemical.Obsolete.Components.Solution solution(temperature_start=310.15) annotation (Placement(transformation(extent={{-100,-56},{100,32}})));

        Chemical.Obsolete.Components.Reaction quaternaryForm(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

        Chemical.Obsolete.Components.Substance O2_free(
          substanceData(DfG=DfG_O2, DfH=-11700),
          use_mass_start=false,
          amountOfSubstance_start=initialO2) annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
        Modelica.Blocks.Sources.ContinuousClock oxygenSource(offset=1000)
          annotation (Placement(transformation(extent={{-78,48},{-58,68}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance oxygen_in_air(
          usePartialPressureInput=true,
          substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
          Temperature=310.15) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-36,58})));
        Chemical.Obsolete.Components.GasSolubility partialPressure1(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,32})));

        Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
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
        Chemical.Obsolete.Sources.ExternalMoleFraction H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          MoleFraction=initialH,
          Temperature=310.15) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={-12,-18})));
        Chemical.Obsolete.Components.Substance CO2_free(
          substanceData(DfG=DfG_CO2, DfH=-412900),
          use_mass_start=false,
          amountOfSubstance_start=initialCO2) annotation (Placement(transformation(extent={{86,-8},{66,12}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO2_gas(
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
          PartialPressure(displayUnit="Pa") = pCO2,
          Temperature=310.15) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={62,60})));
        Chemical.Obsolete.Components.GasSolubility partialPressure2(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,32})));

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
</html>", info="<html>
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
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
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

        Chemical.Obsolete.Components.Speciation speciation(NumberOfSubunits=N) annotation (Placement(transformation(extent={{-18,-72},{2,-52}})));
         // AmountOfSubstance_start=4e-11)
        // AmountOfSubstance_start=totalAmountOfHemoglobin)
        Chemical.Obsolete.Components.Substance OxyHm[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_O2 + RT*log(Ko25) + DfG_selectedForm/N, DfH=DfH_O2 - Ho + DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=(initialO2/Ko37)*initialHb/(1 + initialO2/Ko37 + initialCO/Kco37)) "Oxygenated subunit"
          annotation (Placement(transformation(extent={{-88,14},{-68,34}})));

        Chemical.Obsolete.Components.Reaction o[4](
          each nS=1,
          each nP=2,
          each KC=KC) annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
        Chemical.Obsolete.Components.Substance DeoxyHm[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=initialHb/(1 + initialO2/Ko37 + initialCO/Kco37)) "Deoxygenated subunit"
          annotation (Placement(transformation(extent={{-8,12},{-28,32}})));

        Chemical.Obsolete.Interfaces.SolutionPort solution
          annotation (Placement(transformation(extent={{-50,-82},{-32,-62}}), iconTransformation(extent={{-50,-90},{-30,-70}})));
        Chemical.Obsolete.Interfaces.SubstancePort_b O2
          annotation (Placement(transformation(extent={{-28,32},{-8,52}}), iconTransformation(extent={{-90,70},{-70,90}})));
        Chemical.Obsolete.Interfaces.SubstancePort_a selectedForm
          annotation (Placement(transformation(extent={{26,-82},{46,-62}}), iconTransformation(extent={{30,-90},{50,-70}})));
        Chemical.Obsolete.Components.Substance HmAH[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=RT*log(Kh25) + DfG_selectedForm/N, DfH=-Hh + DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=initialH*initialHb/(Kh37 + initialH)) "Protonated h site of subunit in quaternary structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{10,12},{30,32}})));

        Chemical.Obsolete.Components.Reaction h[4](
          each nS=1,
          each nP=2,
          each KC=KC) annotation (Placement(transformation(extent={{36,32},{56,12}})));
        Chemical.Obsolete.Components.Substance HmA[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=Kh37*initialHb/(Kh37 + initialH)) "Deprotonated h site of subunit in quaternary structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{92,14},{72,34}})));

        Chemical.Obsolete.Components.Substance HmNH3[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=RT*log(Kz25) + DfG_selectedForm/N, DfH=-Hz + DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=(initialH^2)*initialHb/(initialH^2 + initialH*Kz37 + Kz37*Kc37*initialCO2))
          "Protonated z site of subunit in quaternary structure of hemoglobin tetramer" annotation (Placement(transformation(extent={{-84,-42},{-64,-22}})));

        Chemical.Obsolete.Components.Reaction z[4](
          each nS=1,
          each nP=2,
          each KC=KC) annotation (Placement(transformation(extent={{-54,-42},{-34,-22}})));
        Chemical.Obsolete.Components.Substance HmNH2[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_selectedForm/N, DfH=DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=initialH*Kz37*initialHb/(initialH^2 + initialH*Kz37 + Kz37*Kc37*initialCO2))
          "Deprotonated z site of subunit in quaternary structure of hemoglobin tetramer" annotation (Placement(transformation(extent={{12,-44},{-8,-24}})));

        Chemical.Obsolete.Components.Reaction c[4](
          each nP=2,
          each KC=KC,
          each nS=2) annotation (Placement(transformation(extent={{20,-42},{40,-22}})));
        Chemical.Obsolete.Interfaces.SubstancePort_b CO2
          annotation (Placement(transformation(extent={{-8,-26},{12,-6}}), iconTransformation(extent={{10,70},{30,90}})));
        Chemical.Obsolete.Components.Substance HmNHCOO[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfG=DfG_CO2 - RT*log(Kc25) + DfG_selectedForm/N, DfH=DfH_CO2 + Hc + DfH_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=(Kz37*Kc37*initialCO2)*initialHb/(initialH^2 + initialH*Kz37 + Kz37*Kc37*initialCO2))
          "Carboxylated c site of subunit in quaternary structure of hemoglobin tetramer" annotation (Placement(transformation(extent={{70,-44},{50,-24}})));

        Obsolete.Interfaces.SubstancePort_b H annotation (Placement(transformation(extent={{50,-12},{70,8}}), iconTransformation(extent={{70,70},{90,90}})));
        Chemical.Obsolete.Components.Substance COHm[4](
          redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible,
          each substanceData(DfH=DfH_CO - Hco + DfH_selectedForm/N, DfG=DfG_CO + RT*log(Kco25) + DfG_selectedForm/N),
          each use_mass_start=false,
          each amountOfSubstance_start=(initialCO/Kco37)*initialHb/(1 + initialO2/Ko37 + initialCO/Kco37)) "Subunit with Carbon Monoxide"
          annotation (Placement(transformation(extent={{78,44},{58,64}})));

        Chemical.Obsolete.Components.Reaction o1[4](
          each nS=1,
          each nP=2,
          each KC=KC) annotation (Placement(transformation(extent={{46,44},{26,64}})));
        Chemical.Obsolete.Interfaces.SubstancePort_b CO
          annotation (Placement(transformation(extent={{-12,62},{8,82}}), iconTransformation(extent={{-50,70},{-30,90}})));
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
</html>", info="<html>
<p>M. Matej&aacute;k, T. Kulh&aacute;nek, and S. Matou&scaron;ek, &quot;Adair-based hemoglobin equilibrium with oxygen, carbon dioxide and hydrogen ion activity,&quot; Scandinavian Journal of Clinical &amp; Laboratory Investigation, pp. 1-8, 2015.</p>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
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

        Chemical.Obsolete.Components.Solution solution(temperature_start=310.15) annotation (Placement(transformation(extent={{-100,-56},{100,32}})));

        Chemical.Obsolete.Components.Reaction quaternaryForm(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

        Chemical.Obsolete.Components.Substance O2_free(
          substanceData(DfG=DfG_O2, DfH=-11700),
          use_mass_start=false,
          amountOfSubstance_start=initialO2) annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
        Modelica.Blocks.Sources.ContinuousClock oxygenSource(offset=2000)
          annotation (Placement(transformation(extent={{-78,48},{-58,68}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance oxygen_in_air(
          usePartialPressureInput=true,
          substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
          Temperature=310.15) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-36,58})));
        Chemical.Obsolete.Components.GasSolubility partialPressure1(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,32})));

        Obsolete.Components.Substance H2O(substanceData=Obsolete.Substances.Water_liquid(), mass_start=1)
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
        Chemical.Obsolete.Sources.ExternalMoleFraction H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          MoleFraction=initialH,
          Temperature=310.15) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={-12,-18})));
        Chemical.Obsolete.Components.Substance CO2_free(
          substanceData(DfG=DfG_CO2, DfH=-412900),
          use_mass_start=false,
          amountOfSubstance_start=initialCO2) annotation (Placement(transformation(extent={{86,-8},{66,12}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO2_gas(
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
          PartialPressure(displayUnit="kPa") = 5330,
          Temperature=310.15) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={62,60})));
        Chemical.Obsolete.Components.GasSolubility partialPressure2(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,32})));

        Real sCO "Hemoglobin carbon monoxide saturation";
        Real sO2 "Hemoglobin oxygen saturation";
        Real sCO2 "Hemoglobin carbon dioxide saturation";
        Real dH "Hemoglobin charge change caused by binding of Bohr's protons";
        Chemical.Obsolete.Components.Substance CO_free(
          substanceData=Chemical.Obsolete.Substances.CarbonMonoxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=initialCO) annotation (Placement(transformation(extent={{-92,8},{-72,28}})));
        Chemical.Obsolete.Components.GasSolubility partialPressure3(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={26,32})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO_gas(
          substanceData=Chemical.Obsolete.Substances.CarbonMonoxide_gas(),
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
</html>", info="<html>
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
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
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

        Chemical.Obsolete.Components.Solution solution(temperature_start=310.15) annotation (Placement(transformation(extent={{-96,-60},{104,28}})));

        Chemical.Obsolete.Components.Reaction quaternaryForm(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

        Chemical.Obsolete.Components.Substance O2_free(
          substanceData(DfG=DfG_O2, DfH=-11700),
          use_mass_start=false,
          amountOfSubstance_start=initialO2) annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance oxygen_in_air(
          usePartialPressureInput=false,
          substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
          PartialPressure(displayUnit="mmHg") = 11999.01486735,
          Temperature=310.15) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-36,58})));
        Chemical.Obsolete.Components.GasSolubility partialPressure1(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,32})));

        Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
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
        Chemical.Obsolete.Sources.ExternalMoleFraction H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          MoleFraction=initialH,
          useMoleFractionInput=true,
          Temperature=310.15) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              origin={-12,-18},
              rotation=270)));
        Chemical.Obsolete.Components.Substance CO2_free(
          substanceData(DfG=DfG_CO2, DfH=-412900),
          use_mass_start=false,
          amountOfSubstance_start=initialCO2) annotation (Placement(transformation(extent={{86,-8},{66,12}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO2_gas(
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
          PartialPressure(displayUnit="mmHg") = 5332.8954966,
          Temperature=310.15) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={62,60})));
        Chemical.Obsolete.Components.GasSolubility partialPressure2(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,32})));

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
</html>", info="<html>
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
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
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

        Chemical.Obsolete.Components.Solution solution(temperature_start=310.15) annotation (Placement(transformation(extent={{-96,-60},{104,28}})));

        Chemical.Obsolete.Components.Reaction quaternaryForm(
          nS=1,
          nP=1,
          KC=KC) annotation (Placement(transformation(extent={{-22,-52},{-2,-32}})));

        Chemical.Obsolete.Components.Substance O2_free(
          substanceData(DfG=DfG_O2, DfH=-11700),
          use_mass_start=false,
          amountOfSubstance_start=initialO2) annotation (Placement(transformation(extent={{-76,-12},{-56,8}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance oxygen_in_air(
          usePartialPressureInput=false,
          substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
          PartialPressure(displayUnit="mmHg") = 11999.01486735,
          Temperature=310.15) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-36,58})));
        Chemical.Obsolete.Components.GasSolubility partialPressure1(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,32})));

        Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
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
        Chemical.Obsolete.Sources.ExternalMoleFraction H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          MoleFraction=initialH,
          useMoleFractionInput=false,
          Temperature=310.15) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              origin={-12,-18},
              rotation=270)));
        Chemical.Obsolete.Components.Substance CO2_free(
          substanceData(DfG=DfG_CO2, DfH=-412900),
          use_mass_start=false,
          amountOfSubstance_start=initialCO2) annotation (Placement(transformation(extent={{86,-8},{66,12}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance CO2_gas(
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
          usePartialPressureInput=true,
          PartialPressure(displayUnit="kPa") = 5330,
          Temperature=310.15) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={62,60})));
        Chemical.Obsolete.Components.GasSolubility partialPressure2(KC=KC) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,32})));

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
</html>", info="<html>
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
</html>"),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
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

        Chemical.Obsolete.Sensors.DissociationCoefficient dissociationCoefficient(nS=1, nP=1) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Chemical.Obsolete.Sources.PureSubstance A(redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible, substanceData(DfG=0))
          annotation (Placement(transformation(extent={{-56,-10},{-36,10}})));
        Chemical.Obsolete.Sources.PureSubstance B(redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.Incompressible, substanceData(DfG=-R*T_25degC*
                log(K))) annotation (Placement(transformation(extent={{60,-10},{40,10}})));
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
</html>"),experiment(StopTime=0.001));
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

        Chemical.Obsolete.Sources.PureSubstance A annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
        Chemical.Obsolete.Sensors.DissociationCoefficient reaction(nS=2, nP=1) annotation (Placement(transformation(extent={{4,-8},{24,12}})));
        Chemical.Obsolete.Sources.PureSubstance B annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
        Chemical.Obsolete.Sources.PureSubstance C(substanceData(DfG=-R*T_25degC*log(Kx))) annotation (Placement(transformation(extent={{68,-8},{48,12}})));

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
</html>"),experiment(StopTime=0.001));
      end SimpleReaction2;

      model SimpleReaction2_Get_DfG
        "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
         extends Modelica.Icons.Example;

        Chemical.Obsolete.Sources.PureSubstance A annotation (Placement(transformation(extent={{-28,42},{-8,62}})));
        Chemical.Obsolete.Sensors.DissociationCoefficient reaction(nP=1, nS=2) annotation (Placement(transformation(extent={{10,32},{30,52}})));
        Chemical.Obsolete.Sources.PureSubstance B annotation (Placement(transformation(extent={{-28,16},{-8,36}})));

        Modelica.Blocks.Math.InverseBlockConstraints inverseBlockConstraints
          annotation (Placement(transformation(extent={{-42,-80},{82,80}})));
        Chemical.Obsolete.Sources.ExternalElectroChemicalPotential C(usePotentialInput=true) annotation (Placement(transformation(extent={{60,32},{40,52}})));
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
</html>"),experiment(StopTime=0.001));
      end SimpleReaction2_Get_DfG;

      model StandardElectrochemicalCell
        "Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential "
       extends Modelica.Icons.Example;

        Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-90,-40},{-46,68}})));

        Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{60,-40},{96,70}})));

        Chemical.Obsolete.Sources.PureSubstance Ag(substanceData=Chemical.Obsolete.Substances.Silver_solid())
          annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
        Chemical.Obsolete.Sources.PureSubstance Cl(substanceData=Chemical.Obsolete.Substances.Chloride_aqueous())
          annotation (Placement(transformation(extent={{-8,-36},{-28,-16}})));
        Chemical.Obsolete.Sources.PureSubstance AgCl(substanceData=Chemical.Obsolete.Substances.SilverChloride_solid())
          annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance H2(
          substanceData=Chemical.Obsolete.Substances.Hydrogen_gas(),
          PartialPressure=100000,
          TotalPressure=100000) annotation (Placement(transformation(extent={{24,32},{44,52}})));
        Chemical.Obsolete.Sources.PureSubstance H(substanceData=Chemical.Obsolete.Substances.Proton_aqueous())
          annotation (Placement(transformation(extent={{18,-36},{38,-16}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-6,64},{14,84}})));
        Chemical.Obsolete.Components.Reaction electrodeReaction(
          nP=2,
          p={2,2},
          nS=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={52,6})));
        Chemical.Obsolete.Components.Reaction electrodeReaction1(nS=2, nP=2)
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-40,6})));
        Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
        Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{86,-26},{66,-6}})));

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
</html>",   revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end StandardElectrochemicalCell;

      model StandardLeadAcidPotential
        "Standard potential of the lead acid battery"
       extends Modelica.Icons.Example;

        Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{54,-46},{92,62}})));

        Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-94,-50},{-56,58}})));

        Chemical.Obsolete.Sources.PureSubstance Pb(substanceData=Chemical.Obsolete.Substances.Lead_solid())
          annotation (Placement(transformation(extent={{84,-34},{64,-14}})));
        Chemical.Obsolete.Sources.PureSubstance HSO4(substanceData=Chemical.Obsolete.Substances.HydrogenSulfate_aqueous())
          annotation (Placement(transformation(extent={{-22,-58},{-2,-38}})));
        Chemical.Obsolete.Sources.PureSubstance PbSO4_(substanceData=Chemical.Obsolete.Substances.LeadSulfate_solid())
          annotation (Placement(transformation(extent={{84,4},{64,24}})));
        Chemical.Obsolete.Sources.PureSubstance H(substanceData=Chemical.Obsolete.Substances.Proton_aqueous())
          annotation (Placement(transformation(extent={{6,-28},{26,-8}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-6,60},{14,80}})));
        Chemical.Obsolete.Components.Reaction electrodeReaction(
          nP=2,
          nS=4,
          s={1,1,3,2},
          p={1,2}) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-42,14})));
        Chemical.Obsolete.Components.Reaction electrodeReaction1(
          nS=2,
          nP=3,
          p={1,1,2}) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,14})));

        Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{84,32},{64,52}})));
        Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{-86,-12},{-66,8}})));
        Chemical.Obsolete.Sources.PureSubstance PbO2(substanceData=Chemical.Obsolete.Substances.LeadDioxide_solid())
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-74,-30})));
        Chemical.Obsolete.Sources.PureSubstance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid())
          annotation (Placement(transformation(extent={{-2,-10},{-22,10}})));
        Chemical.Obsolete.Sources.PureSubstance PbSO4(substanceData=Chemical.Obsolete.Substances.LeadSulfate_solid())
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-74,32})));

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

        Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-90,-40},{-46,68}})));

        Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{60,-40},{96,70}})));

        Chemical.Obsolete.Sources.PureSubstance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid())
          annotation (Placement(transformation(extent={{-8,-36},{-28,-16}})));
        Chemical.Obsolete.Sources.PureSubstance O2(redeclare package stateOfMatter = Obsolete.Interfaces.IdealGas, substanceData=
              Chemical.Obsolete.Substances.Oxygen_gas()) annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
        Chemical.Obsolete.Sources.ExternalIdealGasSubstance H2(
          substanceData=Chemical.Obsolete.Substances.Hydrogen_gas(),
          PartialPressure=100000,
          TotalPressure=100000) annotation (Placement(transformation(extent={{24,32},{44,52}})));
        Chemical.Obsolete.Sources.PureSubstance H(substanceData=Chemical.Obsolete.Substances.Proton_aqueous())
          annotation (Placement(transformation(extent={{18,-36},{38,-16}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-6,64},{14,84}})));
        Chemical.Obsolete.Components.Reaction electrodeReaction(
          nP=2,
          p={2,2},
          nS=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={52,6})));
        Chemical.Obsolete.Components.Reaction electrodeReaction1(
          s={4},
          p={2,8,8},
          nS=1,
          nP=3) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-40,6})));
        Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
        Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{86,-26},{66,-6}})));

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
</html>",   revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end StandardElectrochemicalCell2;
    end CheckSubstancesData;

    model WaterElectrolysis "Water electrolysis"
      extends Modelica.Icons.Example;
      Chemical.Obsolete.Components.Substance O2_gas(
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=0.001) annotation (Placement(transformation(extent={{-22,-6},{-2,14}})));

      Chemical.Obsolete.Components.Substance H2_gas(
        redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.IdealGas,
        substanceData=Chemical.Obsolete.Substances.Hydrogen_gas(),
        use_mass_start=false,
        amountOfSubstance_start=0.001) annotation (Placement(transformation(extent={{16,-6},{36,14}})));
      Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{58,-78},{92,30}})));
      Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-90,-80},{-56,28}})));
      Chemical.Obsolete.Components.Solution water(temperature_start=310.15) annotation (Placement(transformation(extent={{-28,-80},{18,-46}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-42,70},{-22,90}})));
      Chemical.Obsolete.Components.Reaction reaction(
        nS=2,
        s={2,4},
        nP=3,
        p={2,1,4}) annotation (Placement(transformation(
            extent={{11,11},{-11,-11}},
            rotation=180,
            origin={-23,-31})));
      Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{84,-38},{64,-18}})));
      Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{-84,-34},{-64,-14}})));
    Modelica.Electrical.Analog.Basic.Resistor resistor(R=1)
      annotation (Placement(transformation(extent={{-36,38},{-16,58}})));
    Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
      annotation (Placement(transformation(extent={{-66,38},{-46,58}})));
      Chemical.Obsolete.Components.Solution air(redeclare package stateOfMatter = Chemical.Obsolete.Interfaces.IdealGas)
        annotation (Placement(transformation(extent={{-40,-16},{50,26}})));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
        annotation (Placement(transformation(extent={{18,38},{-2,58}})));
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{38,26},{58,46}})));
      Obsolete.Components.Substance liquidWater(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
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

      Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
      Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{62,-50},{96,50}})));

      Chemical.Obsolete.Components.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-30,-60},{38,6}})));

      Chemical.Obsolete.Sources.ExternalIdealGasSubstance H2(
        substanceData=Chemical.Obsolete.Substances.Hydrogen_gas(),
        PartialPressure=100000,
        TotalPressure=100000) annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Chemical.Obsolete.Components.Substance H(
        substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-6,-22},{14,-2}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Chemical.Obsolete.Components.Reaction electrodeReaction(
        p={2,2},
        nS=1,
        nP=2) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={52,6})));
      Chemical.Obsolete.Components.Reaction electrodeReaction1(
        s={4},
        p={2,8,8},
        nS=1,
        nP=3) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-40,0})));

      Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                                 //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                  //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{62,54},{82,74}})));
      Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{-6,-54},{14,-34}})));
      Obsolete.Sources.ExternalIdealGasSubstance O2_(
        substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
        PartialPressure=100000,
        TotalPressure=100000) annotation (Placement(transformation(extent={{-6,32},{-26,52}})));
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
        Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
        Obsolete.Components.Substance CH3COOH(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.AceticAcid_aqueous(),
          mass_start=0.001) "Acetic acid" annotation (Placement(transformation(extent={{-72,30},{-52,50}})));
        Obsolete.Components.Substance CH4(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.Methan_aqueous(),
          mass_start=0.001) "Methan" annotation (Placement(transformation(extent={{26,48},{46,68}})));
        Obsolete.Components.Substance CO2(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
          mass_start=0.001) annotation (Placement(transformation(extent={{56,10},{76,30}})));
        Obsolete.Components.Substance Water(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.Water_liquid(),
          mass_start=1) annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));
        Obsolete.Components.Reaction reaction(nS=1, nP=2) "Acetoclastic (heterotrophic) methanogenesis"
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
        Obsolete.Components.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
        Obsolete.Components.Substance CH4(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.Methan_aqueous(),
          mass_start=0.001) "Methan" annotation (Placement(transformation(extent={{32,70},{52,90}})));
        Obsolete.Components.Substance CO2(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
          mass_start=0.001) annotation (Placement(transformation(extent={{-68,40},{-48,60}})));
        Obsolete.Components.Substance H2O(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.Water_liquid(),
          mass_start=1) annotation (Placement(transformation(extent={{36,36},{56,56}})));
        Obsolete.Components.Reaction reaction(
          KC=1e-7,
          s={4,1},
          p={1,2},
          nS=2,
          nP=2) "Hydrogenotrophic (autotrophic) methanogenesis" annotation (Placement(transformation(extent={{-18,50},{2,70}})));
        Obsolete.Components.Substance H2(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.Hydrogen_aqueous(),
          mass_start=0.001) annotation (Placement(transformation(extent={{-74,66},{-54,86}})));
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

        Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
        Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{62,-50},{96,50}})));

        Chemical.Obsolete.Components.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-30,-96},{38,6}})));

        Chemical.Obsolete.Sources.ExternalIdealGasSubstance CH4(
          substanceData=Chemical.Obsolete.Substances.Methan_gas(),
          PartialPressure=100000,
          TotalPressure=100000) annotation (Placement(transformation(extent={{22,26},{42,46}})));
        Chemical.Obsolete.Components.Substance H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-6,-22},{14,-2}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-20,62},{0,82}})));
        Chemical.Obsolete.Components.Reaction electrodeReaction(
          s={1,8,8},
          p={1,2},
          nP=2,
          nS=3) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={52,6})));
        Chemical.Obsolete.Components.Reaction electrodeReaction1(
          s={4},
          p={2,8,8},
          nS=1,
          nP=3) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-40,0})));

        Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                                   //(substanceData=Chemical.Examples.Substances.Electrone_solid())
        Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                    //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{46,52},{66,72}})));
        Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-6,-80},{14,-60}})));
        Obsolete.Sources.ExternalIdealGasSubstance O2_(
          substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
          PartialPressure=100000,
          TotalPressure=100000) annotation (Placement(transformation(extent={{-6,32},{-26,52}})));
        Obsolete.Sources.ExternalIdealGasSubstance CO2(
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_gas(),
          PartialPressure=100000,
          TotalPressure=100000) annotation (Placement(transformation(extent={{88,-90},{68,-70}})));
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

        Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-86,26},{-46,72}})));
        Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{58,26},{94,72}})));

        Chemical.Obsolete.Components.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-90,-90},{92,14}})));

        Chemical.Obsolete.Components.Substance H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-84,-14},{-64,6}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-6,64},{14,84}})));
        Chemical.Obsolete.Components.Reaction electrodeReaction(
          p={2,2},
          nP=2,
          nS=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={50,16})));
        Chemical.Obsolete.Components.Reaction electrodeReaction1(
          s={4},
          p={2,8,8},
          nS=1,
          nP=3) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-36,28})));

        Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
                                   //(substanceData=Chemical.Examples.Substances.Electrone_solid())
        Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{88,38},{68,58}})));
                                    //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{22,54},{42,74}})));
        Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{78,-28},{58,-8}})));
        Obsolete.Sources.ExternalIdealGasSubstance O2_(
          substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
          PartialPressure=100000,
          TotalPressure=100000) annotation (Placement(transformation(extent={{0,36},{-20,56}})));
        Obsolete.Components.Substance CH4(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.Methan_aqueous(),
          mass_start=0.001) "Methan" annotation (Placement(transformation(extent={{74,-64},{54,-44}})));
        Obsolete.Components.Substance CO2(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
          mass_start=0.1) annotation (Placement(transformation(extent={{-72,-76},{-52,-56}})));
        Obsolete.Components.Reaction reaction(
          KC=1e-7,
          s={4,1},
          p={1,2},
          nS=2,
          nP=2) "Hydrogenotrophic (autotrophic) methanogenesis" annotation (Placement(transformation(extent={{0,-54},{20,-34}})));
        Obsolete.Components.Substance H2(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.Hydrogen_aqueous(),
          mass_start=5e-11) annotation (Placement(transformation(extent={{-76,-44},{-56,-24}})));
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

        Chemical.Obsolete.Components.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-86,26},{-46,72}})));
        Chemical.Obsolete.Components.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{58,26},{94,72}})));

        Chemical.Obsolete.Components.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-90,-90},{92,14}})));

        Chemical.Obsolete.Components.Substance H(
          substanceData=Chemical.Obsolete.Substances.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-4) annotation (Placement(transformation(extent={{-84,-14},{-64,6}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-6,64},{14,84}})));
        Chemical.Obsolete.Components.Reaction electrodeReaction1(
          s={4},
          p={2,8,8},
          nS=1,
          nP=3) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-36,28})));

        Chemical.Obsolete.Components.ElectronTransfer electrone annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
                                   //(substanceData=Chemical.Examples.Substances.Electrone_solid())
        Chemical.Obsolete.Components.ElectronTransfer electrone1 annotation (Placement(transformation(extent={{88,38},{68,58}})));
                                    //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{22,54},{42,74}})));
        Obsolete.Components.Substance H2O(substanceData=Chemical.Obsolete.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{78,-28},{58,-8}})));
        Obsolete.Sources.ExternalIdealGasSubstance O2_(
          substanceData=Chemical.Obsolete.Substances.Oxygen_gas(),
          PartialPressure=100000,
          TotalPressure=100000) annotation (Placement(transformation(extent={{0,36},{-20,56}})));
        Obsolete.Components.Substance AcAc(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.AceticAcid_aqueous(),
          mass_start=0.001) "Acetic Acid" annotation (Placement(transformation(extent={{60,-60},{40,-40}})));
        Obsolete.Components.Substance CO2(
          redeclare package stateOfMatter = Obsolete.Interfaces.Incompressible,
          substanceData=Chemical.Obsolete.Substances.CarbonDioxide_aqueous(),
          mass_start=1) annotation (Placement(transformation(extent={{-10,-34},{10,-14}})));
        Obsolete.Components.Reaction reaction(
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
end Obsolete;
