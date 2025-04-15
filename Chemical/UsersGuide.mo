within Chemical;
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
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Solution\">Chemical solution</a></p><p>The solution is the base component of each model, because it defines the conditions of the electro-chemical processes. It integrates the total amount of substance (called amount of solution), heat, charge, entropy, volume and others from each substances to present the base properties such as temperature, pressure, electric potential and others. The usage is very simple - just connect each chemical substance with its chemical solution using their <a href=\"modelica://Chemical.Interfaces.SolutionPort\">SolutionPort</a>.</p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Substance1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.Substance\">Chemical substance</a></p><p>The chemical substance integrates the amount of the chemical substance and from the properties of the connected solution it presents the electro-chemical potential of the substance using the <a href=\"modelica://Chemical.Obsolete.Interfaces.SubstancePort\">SubstancePort</a>.</p><p>There are two basic <a href=\"modelica://Chemical.Interfaces.StateOfMatter\">states of matter</a>: <a href=\"modelica://Chemical.Interfaces.IdealGas\">ideal gas</a> and <a href=\"modelica://Chemical.Interfaces.Incompressible\">incompressible</a> substance. However, the user can easily (re)define their own state of matter by inserting the correct expressions for the pure substance <a href=\"modelica://Chemical.Interfaces.StateOfMatter.activityCoefficient\">activity coefficient</a>, <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarVolumePure\">molar volume</a>, <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarEntropyPure\">molar entropy</a> and <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarEnthalpyElectroneutral\">molar enthalpy</a>, based on the current solution state (temperature, pressure, electric potential and ionic strength) and the <a href=\"modelica://Chemical.Interfaces.StateOfMatter.SubstanceData\">substance data</a>. The object-oriented design allows users to define the substance data record as part of the state of matter package. Users can select substance parameters according to the state of matter, redefining the getter functions of substance properties.</p><p>The examples work with ideal gases in case of all gaseous substance and incompressible state of matter in case of liquid or solid. The definition data are the molar mass of the substance, the number of charges of the substance, the molar heat capacity of the substance at a constant pressure, free formation enthalpy, free formation Gibbs energy and density (if incompressible) &mdash; all at a temperature of 25&deg;C and pressure 1 bar. Since these parameters are usually recorded in chemical tables at this standard conditions. In this manner, more than 35 real chemical <a href=\"modelica://Examples.Substances\">substances</a> in the example package of this chemical library have already been defined. The usage of these predefined substances&rsquo; data is very simple. In the parameter dialog of the chemical substance, the correct record with this data can be selected, as shown in Figure 1.</p><p>This setting is typically the most important setting of each chemical model. All equilibrium coefficients, standard voltages, dissolution coefficients, saturated vapor pressures and so on, are automatically solved using these substance data. As a result, for example, the chemical reaction component only needs to define the stoichiometry coefficients, and the connected substances reach equilibrium at the correct equilibrium coefficient.</p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Reaction1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.Reaction\">Chemical reaction</a></p><p>The chemical reaction component is very general. The dissociation constant of the equilibrium is calculated from substance properties at usual in thermodynamics, for example as definition of <a href=\"http://goldbook.iupac.org/S05915.html\">UIPAC</a>. For example if we want to define <a href=\"modelica://Examples.SimpleReaction\">simple reaction A&lt;-&gt;B</a> with dissociation constant [B]/[A]=2 then it must be the difference between Gibbs energies of formation equal to B.DfG - A.DfG = - R * T * ln(2). Without lost of generality it is possible to select some substances as reference and give them the zero Gibbs energy of formation. The next substances created by some chemical process can be expressed from them such as example of <a href=\"modelica://Examples.Hemoglobin.Allosteric_Hemoglobin_MWC\">alosteric hemoglobin</a> calculation. The kinetics of the chemical reaction is different as usual. However the most of processes can be recalculated with sufficient precision, for example the <a href=\"modelica://Examples.EnzymeKinetics\">Michaelic-Menton</a> can be recalculated with precision of 1.5% of maximal rate. </p></td>
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
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Obsolete.Components.Speciation\">Chemical speciation</a></p><p>The chemical speciation is for macromolecule composed with independent subunits is specific conformations. For example the hemoglobin is tetramer, which can be in two conformation: relaxed and tensed. In each of this conformation it has different afinities (different dissociation constant) for binding oxygen in each of four independent subunits. This alosteric effect can be modeled using speciation such as in <a href=\"modelica://Examples.Hemoglobin.Allosteric_Hemoglobin2_MWC\">Allosteric_Hemoglobin2_MWC</a>. However the result should be the same as using the detailed reaction model <a href=\"modelica://Examples.Hemoglobin.Allosteric_Hemoglobin_MWC\">Allosteric_Hemoglobin_MWC</a>.</p></td>
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
<td valign=\"middle\"><p><br><a href=\"modelica://Chemical.Interfaces.SolutionPort\">Chemical.Interfaces.SolutionPort</a></p></td>
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
