within ;
package Chemical "Library of Electro-Chemical models (chemical reactions, diffusions, membrane channels, gas dissolutions, electrochemical cells, ...)"
  package UsersGuide "User's Guide"
    extends Modelica.Icons.Information;

  class Overview "Overview"
    extends Modelica.Icons.Information;

   annotation (Documentation(info="<html>
<p>The Chemical library can describe the following phenomena.</p>
<table cellspacing=\"0\" cellpadding=\"2\" border=\"1\"><tr>
<td><p align=\"center\"><h4>Chemical Components</h4></p></td>
<td><p align=\"center\"><h4>Description</h4></p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Solution1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Solution\">Chemical solution</a></p><p>The solution is the base component of each model, because it defines the conditions of the electro-chemical processes. It integrates the total amount of substance (called amount of solution), heat, charge, entropy, volume and others from each substances to present the base properties such as temperature, pressure, electric potential and others. The usage is very simple - just connect each chemical substance with its chemical solution using their <a href=\"modelica://Chemical.Interfaces.SolutionPort\">SolutionPort</a>.</p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Substance1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Substance\">Chemical substance</a></p><p>The chemical substance integrates the amount of the chemical substance and from the properties of the connected solution it presents the electro-chemical potential of the substance using the <a href=\"modelica://Chemical.Interfaces.ChemicalPort\">SubstancePort</a>.</p><p>There are two basic <a href=\"modelica://Chemical.Interfaces.StateOfMatter\">states of matter</a>: <a href=\"modelica://Chemical.Interfaces.IdealGas\">ideal gas</a> and <a href=\"modelica://Chemical.Interfaces.Incompressible\">incompressible</a> substance. However, the user can easily (re)define their own state of matter by inserting the correct expressions for the pure substance <a href=\"modelica://Chemical.Interfaces.StateOfMatter.activityCoefficient\">activity coefficient</a>, <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarVolumePure\">molar volume</a>, <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarEntropyPure\">molar entropy</a> and <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarEnthalpyElectroneutral\">molar enthalpy</a>, based on the current solution state (temperature, pressure, electric potential and ionic strength) and the <a href=\"modelica://Chemical.Interfaces.StateOfMatter.SubstanceData\">substance data</a>. The object-oriented design allows users to define the substance data record as part of the state of matter package. Users can select substance parameters according to the state of matter, redefining the getter functions of substance properties.</p><p>The examples work with ideal gases in case of all gaseous substance and incompressible state of matter in case of liquid or solid. The definition data are the molar mass of the substance, the number of charges of the substance, the molar heat capacity of the substance at a constant pressure, free formation enthalpy, free formation Gibbs energy and density (if incompressible) &mdash; all at a temperature of 25&deg;C and pressure 1 bar. Since these parameters are usually recorded in chemical tables at this standard conditions. In this manner, more than 35 real chemical <a href=\"modelica://Chemical.Examples.Substances\">substances</a> in the example package of this chemical library have already been defined. The usage of these predefined substances&rsquo; data is very simple. In the parameter dialog of the chemical substance, the correct record with this data can be selected, as shown in Figure 1.</p><p>This setting is typically the most important setting of each chemical model. All equilibrium coefficients, standard voltages, dissolution coefficients, saturated vapor pressures and so on, are automatically solved using these substance data. As a result, for example, the chemical reaction component only needs to define the stoichiometry coefficients, and the connected substances reach equilibrium at the correct equilibrium coefficient.</p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Reaction1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Reaction\">Chemical reaction</a></p><p>The chemical reaction component is very general. The dissociation constant of the equilibrium is calculated from substance properties at usual in thermodynamics, for example as definition of <a href=\"http://goldbook.iupac.org/S05915.html\">UIPAC</a>. For example if we want to define <a href=\"modelica://Chemical.Examples.SimpleReaction\">simple reaction A&LT;-&GT;B</a> with dissociation constant [B]/[A]=2 then it must be the difference between Gibbs energies of formation equal to B.DfG - A.DfG = - R * T * ln(2). Without lost of generality it is possible to select some substances as reference and give them the zero Gibbs energy of formation. The next substances created by some chemical process can be expressed from them such as example of <a href=\"modelica://Chemical.Examples.Hemoglobin.Allosteric_Hemoglobin_MWC\">alosteric hemoglobin</a> calculation. The kinetics of the chemical reaction is different as usual. However the most of processes can be recalculated with sufficient precision, for example the <a href=\"Chemical.Examples.MichaelisMenten\">Michaelic-Menton</a> can be recalculated with precision of 1.5&percnt; of maximal rate. </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Diffusion1.png\"/></p></td>
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Diffusion\">Diffusion</a></p><p>Diffusion is a dynamic chemical process, wich is also equilibrating of electro-chemical potential of the substance. Analogically as in chemical reaction the speed of diffucion can be calculated as coefficient C multiplied by electro-chemical gratient. C can be a parammeter or input expressed from distance, substance and solution properties. </p></td>
</tr>
<tr>
<td valign=\"top\"><p align=\"center\"><img src=\"modelica://Chemical/Resources/Images/UsersGuide/Gassolubility1.png\"/></p></td>
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
<td valign=\"top\"><p><br><h4>potential</h4></p><p>variables</p></td>
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
<td valign=\"middle\"><p><br><a href=\"Chemical.Interfaces.ChemicalPort\">Chemical.Interfaces.ChemicalPort</a> </p><p>ChemicalDefinitionPort, ChemicalUsePort</p></td>
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
</table>
</html>"));
  end Connectors;

    class ModelicaLicense2 "Modelica License 2"

      annotation (Documentation(info="<html>
<p>All files in this directory (Chemical) and in all subdirectories, especially all files that build package &QUOT;Chemical&QUOT; are licensed by <u><b>Marek Matej&aacute;k</b></u> under the <u><b>Modelica License 2 </b></u>(with exception of files &QUOT;Resources/*&QUOT;). </p>
<h4>Licensor:</h4>
<p>Marek Matej&aacute;k,</p>
<p>Hviezdoslavova 632/41,</p>
<p>916 01 Star&aacute; Tur&aacute;, </p>
<p>Slovak Republic, </p>
<p>Charles University in Prague, Czech Republic</p>
<p><br>email: marek@matfyz.cz</p>
<h4>Copyright notices of the files:</h4>
<p>Copyright &copy; 2008-2015, Marek Matejak, Charles University in Prague, First Faculty of Medicine, Institute of Pathological Physiology</p>
<p><br><br>This package with all of its subpackages is released under the &ldquo;Modelica License 2&rdquo; (if not explicitly noted otherwise). </p>
<p><br><a href=\"#The_Modelica_License_2-outline\">The Modelica License 2</a></p>
<p><br><a href=\"#How_to_Apply_the_Modelica_License_2-outline\">How to Apply the Modelica License 2</a></p>
<p><br><a href=\"#Frequently_Asked_Questions-outline\">Frequently Asked Questions</a></p>
<p><br><b></font><font style=\"color: #008000; \">The Modelica License 2</font></b> </p>
<p><b><font style=\"font-size: 10pt; \">Preamble.</b> The goal of this license is that Modelica related model libraries, software, images, documents, data files etc. can be used freely in the original or a modified form, in open source and in commercial environments (as long as the license conditions below are fulfilled, in particular sections 2c) and 2d). The Original Work is provided free of charge and the use is completely at your own risk. Developers of free Modelica packages are encouraged to utilize this license for their work. </p>
<p>The Modelica License applies to any Original Work that contains the following licensing notice adjacent to the copyright notice(s) for this Original Work: </p>
<p><b>Licensed by Marek Matejak under the Modelica License 2</b> </p>
<h4>1. Definitions.</h4>
<p>&ldquo;License&rdquo; is this Modelica License. </p>
<p>&ldquo;Original Work&rdquo; is any work of authorship, including software, images, documents, data files, that contains the above licensing notice or that is packed together with a licensing notice referencing it. </p>
<p>&ldquo;Licensor&rdquo; is the provider of the Original Work who has placed this licensing notice adjacent to the copyright notice(s) for the Original Work. The Original Work is either directly provided by the owner of the Original Work, or by a licensee of the owner. </p>
<p>&ldquo;Derivative Work&rdquo; is any modification of the Original Work which represents, as a whole, an original work of authorship. For the matter of clarity and as examples: </p>
<p>Derivative Work shall not include work that remains separable from the Original Work, as well as merely extracting a part of the Original Work without modifying it. </p>
<p>Derivative Work shall not include (a) fixing of errors and/or (b) adding vendor specific Modelica annotations and/or (c) using a subset of the classes of a Modelica package, and/or (d) using a different representation, e.g., a binary representation. </p>
<p>Derivative Work shall include classes that are copied from the Original Work where declarations, equations or the documentation are modified. </p>
<p>Derivative Work shall include executables to simulate the models that are generated by a Modelica translator based on the Original Work (of a Modelica package). </p>
<p>&ldquo;Modified Work&rdquo; is any modification of the Original Work with the following exceptions: (a) fixing of errors and/or (b) adding vendor specific Modelica annotations and/or (c) using a subset of the classes of a Modelica package, and/or (d) using a different representation, e.g., a binary representation. </p>
<p>&QUOT;Source Code&QUOT; means the preferred form of the Original Work for making modifications to it and all available documentation describing how to modify the Original Work. </p>
<p>&ldquo;You&rdquo; means an individual or a legal entity exercising rights under, and complying with all of the terms of, this License. </p>
<p>&ldquo;Modelica package&rdquo; means any Modelica library that is defined with the &ldquo;<b>package</b></font><font style=\"font-size: 9pt; \">&nbsp;&LT;Name&GT;&nbsp;...&nbsp;</font><font style=\"font-size: 10pt; \">end</font><font style=\"font-size: 9pt; \">&nbsp;&LT;Name&GT;;</font><font style=\"font-size: 10pt; \">&rdquo; Modelica language element. </p>
<p><b>2. Grant of Copyright License.</b> Licensor grants You a worldwide, royalty-free, non-exclusive, sublicensable license, for the duration of the copyright, to do the following: </p>
<p>To reproduce the Original Work in copies, either alone or as part of a collection. </p>
<p>To create Derivative Works according to Section 1d) of this License. </p>
<p>To distribute or communicate to the public copies of the <u>Original Work</u> or a <u>Derivative Work</u> under <u>this License</u>. No fee, neither as a copyright-license fee, nor as a selling fee for the copy as such may be charged under this License. Furthermore, a verbatim copy of this License must be included in any copy of the Original Work or a Derivative Work under this License.</p>
<p>For the matter of clarity, it is permitted A) to distribute or communicate such copies as part of a (possible commercial) collection where other parts are provided under different licenses and a license fee is charged for the other parts only and B) to charge for mere printing and shipping costs. </p>
<p>To distribute or communicate to the public copies of a <u>Derivative Work</u>, alternatively to Section 2c), under <u>any other license</u> of your choice, especially also under a license for commercial/proprietary software, as long as You comply with Sections 3, 4 and 8 below. </p>
<p>For the matter of clarity, no restrictions regarding fees, either as to a copyright-license fee or as to a selling fee for the copy as such apply. </p>
<p>To perform the Original Work publicly. </p>
<p>To display the Original Work publicly. </p>
<p><b>3. Acceptance.</b> Any use of the Original Work or a Derivative Work, or any action according to either Section 2a) to 2f) above constitutes Your acceptance of this License. </p>
<p><b>4. Designation of Derivative Works and of Modified Works. </b>The identifying designation of Derivative Work and of Modified Work must be different to the corresponding identifying designation of the Original Work. This means especially that the (root-level) name of a Modelica package under this license must be changed if the package is modified (besides fixing of errors, adding vendor specific Modelica annotations, using a subset of the classes of a Modelica package, or using another representation, e.g. a binary representation). </p>
<p><b>5. Grant of Patent License.</b> Licensor grants You a worldwide, royalty-free, non-exclusive, sublicensable license, under patent claims owned by the Licensor or licensed to the Licensor by the owners of the Original Work that are embodied in the Original Work as furnished by the Licensor, for the duration of the patents, to make, use, sell, offer for sale, have made, and import the Original Work and Derivative Works under the conditions as given in Section 2. For the matter of clarity, the license regarding Derivative Works covers patent claims to the extent as they are embodied in the Original Work only. </p>
<p><b>6. Provision of Source Code.</b> Licensor agrees to provide You with a copy of the Source Code of the Original Work but reserves the right to decide freely on the manner of how the Original Work is provided.</p>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For the matter of clarity, Licensor might provide only a binary representation of the Original Work. In that case, You may (a) either reproduce the Source Code from the binary representation if this is possible (e.g., by performing a copy of an encrypted Modelica package, if encryption allows the copy operation) or (b) request the Source Code from the Licensor who will provide it to You. </p>
<p><b>7. Exclusions from License Grant.</b> Neither the names of Licensor, nor the names of any contributors to the Original Work, nor any of their trademarks or service marks, may be used to endorse or promote products derived from this Original Work without express prior permission of the Licensor. Except as otherwise expressly stated in this License and in particular in Sections 2 and 5, nothing in this License grants any license to Licensor&rsquo;s trademarks, copyrights, patents, trade secrets or any other intellectual property, and no patent license is granted to make, use, sell, offer for sale, have made, or import embodiments of any patent claims.</p>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;No license is granted to the trademarks of Licensor even if such trademarks are included in the Original Work, except as expressly stated in this License. Nothing in this License shall be interpreted to prohibit Licensor from licensing under terms different from this License any Original Work that Licensor otherwise would have a right to license. </p>
<p><b>8. Attribution Rights.</b> You must retain in the Source Code of the Original Work and of any Derivative Works that You create, all author, copyright, patent, or trademark notices, as well as any descriptive text identified therein as an &QUOT;Attribution Notice&QUOT;. The same applies to the licensing notice of this License in the Original Work. For the matter of clarity, &ldquo;author notice&rdquo; means the notice that identifies the original author(s). </p>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;You must cause the Source Code for any Derivative Works that You create to carry a prominent Attribution Notice reasonably calculated to inform recipients that You have modified the Original Work. </p>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In case the Original Work or Derivative Work is not provided in Source Code, the Attribution Notices shall be appropriately displayed, e.g., in the documentation of the Derivative Work. </p>
<h4>9. Disclaimer of Warranty. </h4>
<p><u>The Original Work is provided under this License on an &QUOT;as is&QUOT; basis and without warranty, either express or implied, including, without limitation, the warranties of non-infringement, merchantability or fitness for a particular purpose. The entire risk as to the quality of the Original Work is with You.</u> This disclaimer of warranty constitutes an essential part of this License. No license to the Original Work is granted by this License except under this disclaimer. </p>
<p><b>10. Limitation of Liability.</b> Under no circumstances and under no legal theory, whether in tort (including negligence), contract, or otherwise, shall the Licensor, the owner or a licensee of the Original Work be liable to anyone for any direct, indirect, general, special, incidental, or consequential damages of any character arising as a result of this License or the use of the Original Work including, without limitation, damages for loss of goodwill, work stoppage, computer failure or malfunction, or any and all other commercial damages or losses. This limitation of liability shall not apply to the extent applicable law prohibits such limitation. </p>
<p><b>11. Termination.</b> This License conditions your rights to undertake the activities listed in Section 2 and 5, including your right to create Derivative Works based upon the Original Work, and doing so without observing these terms and conditions is prohibited by copyright law and international treaty. Nothing in this License is intended to affect copyright exceptions and limitations. This License shall terminate immediately and You may no longer exercise any of the rights granted to You by this License upon your failure to observe the conditions of this license. </p>
<p><b>12. Termination for Patent Action.</b> This License shall terminate automatically and You may no longer exercise any of the rights granted to You by this License as of the date You commence an action, including a cross-claim or counterclaim, against Licensor, any owners of the Original Work or any licensee alleging that the Original Work infringes a patent. This termination provision shall not apply for an action alleging patent infringement through combinations of the Original Work under combination with other software or hardware. </p>
<p><b>13. Jurisdiction.</b> Any action or suit relating to this License may be brought only in the courts of a jurisdiction wherein the Licensor resides and under the laws of that jurisdiction excluding its conflict-of-law provisions. The application of the United Nations Convention on Contracts for the International Sale of Goods is expressly excluded. Any use of the Original Work outside the scope of this License or after its termination shall be subject to the requirements and penalties of copyright or patent law in the appropriate jurisdiction. This section shall survive the termination of this License. </p>
<p><b>14. Attorneys&rsquo; Fees.</b> In any action to enforce the terms of this License or seeking damages relating thereto, the prevailing party shall be entitled to recover its costs and expenses, including, without limitation, reasonable attorneys&apos; fees and costs incurred in connection with such action, including any appeal of such action. This section shall survive the termination of this License. </p>
<p><b>15. Miscellaneous.</b> </p>
<p>If any provision of this License is held to be unenforceable, such provision shall be reformed only to the extent necessary to make it enforceable. </p>
<p>No verbal ancillary agreements have been made. Changes and additions to this License must appear in writing to be valid. This also applies to changing the clause pertaining to written form. </p>
<p>You may use the Original Work in all ways not otherwise restricted or conditioned by this License or by law, and Licensor promises not to interfere with or be responsible for such uses by You. </p>
<p><br><b></font><font style=\"color: #008000; \">How to Apply the Modelica License 2</font></b> </p>
<p><font style=\"font-size: 10pt; \">At the top level of your Modelica package and at every important subpackage, add the following notices in the info layer of the package: </p>
<p>Licensed by &LT;Licensor&GT; under the Modelica License 2</p>
<p>Copyright &copy; &LT;year1&GT;-&LT;year2&GT;, &LT;name of copyright holder(s)&GT;. </p>
<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">http://www.modelica.org/licenses/ModelicaLicense2</a>.</i> </p>
<p>Include a copy of the Modelica License 2 under <b>&LT;library&GT;.UsersGuide.ModelicaLicense2</b> (use <a href=\"http://www.modelica.org/licenses/ModelicaLicense2.mo\">http://www.modelica.org/licenses/ModelicaLicense2.mo</a>). Furthermore, add the list of authors and contributors under <b>&LT;library&GT;.UsersGuide.Contributors</b> or <b>&LT;library&GT;.UsersGuide.Contact</b>. </p>
<p>For example, sublibrary Modelica.Blocks of the Modelica Standard Library may have the following notices: </p>
<p>Licensed by Modelica Association under the Modelica License 2</p>
<p>Copyright &copy; 1998-2008, Modelica Association. </p>
<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">http://www.modelica.org/licenses/ModelicaLicense2</a>.</i> </p>
<p>For C-source code and documents, add similar notices in the corresponding file. </p>
<p>For images, add a &ldquo;readme.txt&rdquo; file to the directories where the images are stored and include a similar notice in this file. </p>
<p>In these cases, save a copy of the Modelica License 2 in one directory of the distribution, e.g., <a href=\"http://www.modelica.org/licenses/ModelicaLicense2.html\">http://www.modelica.org/licenses/ModelicaLicense2.html</a> in directory <b>&LT;library&GT;/Resources/Documentation/ModelicaLicense2.html</b>. </p>
<p><br><b><font style=\"font-size: 6pt; color: #008000; \">Frequently Asked Questions</font></b></p>
<p><font style=\"font-size: 10pt; \">This section contains questions/answer to users and/or distributors of Modelica packages and/or documents under Modelica License 2. Note, the answers to the questions below are not a legal interpretation of the Modelica License 2. In case of a conflict, the language of the license shall prevail. </p>
<p><b></font><font style=\"color: #008000; \">Using or Distributing a Modelica <u>Package</u> under the Modelica License 2</font></b> </p>
<p><b><font style=\"font-size: 10pt; \">What are the main differences to the previous version of the Modelica License?</b></p>
<ol>
<li>Modelica License 1 is unclear whether the licensed Modelica package can be distributed under a different license. Version 2 explicitly allows that &ldquo;Derivative Work&rdquo; can be distributed under any license of Your choice, see examples in Section 1d) as to what qualifies as Derivative Work (so, version 2 is clearer). </li>
<li>If You modify a Modelica package under Modelica License 2 (besides fixing of errors, adding vendor specific Modelica annotations, using a subset of the classes of a Modelica package, or using another representation, e.g., a binary representation), you must rename the root-level name of the package for your distribution. In version 1 you could keep the name (so, version 2 is more restrictive). The reason of this restriction is to reduce the risk that Modelica packages are available that have identical names, but different functionality. </li>
<li>Modelica License 1 states that &ldquo;It is not allowed to charge a fee for the original version or a modified version of the software, besides a reasonable fee for distribution and support&rdquo;. Version 2 has a similar intention for all Original Work under <u>Modelica License 2</u> (to remain free of charge and open source) but states this more clearly as &ldquo;No fee, neither as a copyright-license fee, nor as a selling fee for the copy as such may be charged&rdquo;. Contrary to version 1, Modelica License 2 has no restrictions on fees for Derivative Work that is provided under a different license (so, version 2 is clearer and has fewer restrictions). </li>
<li>Modelica License 2 introduces several useful provisions for the licensee (articles 5, 6, 12), and for the licensor (articles 7, 12, 13, 14) that have no counter part in version 1. </li>
<li>Modelica License 2 can be applied to all type of work, including documents, images and data files, contrary to version 1 that was dedicated for software only (so, version 2 is more general). </li>
</ol>
<h4>Can I distribute a Modelica package (under Modelica License 2) as part of my commercial Modelica modeling and simulation environment?</h4>
<p>Yes, according to Section 2c). However, you are not allowed to charge a fee for this part of your environment. Of course, you can charge for your part of the environment. </p>
<h4>Can I distribute a Modelica package (under Modelica License 2) under a different license?</h4>
<p>No. The license of an unmodified Modelica package cannot be changed according to Sections 2c) and 2d). This means that you cannot <u>sell</u> copies of it, any distribution has to be free of charge. </p>
<h4>Can I distribute a Modelica package (under Modelica License 2) under a different license when I first encrypt the package?</h4>
<p>No. Merely encrypting a package does not qualify for Derivative Work and therefore the encrypted package has to stay under Modelica License 2. </p>
<h4>Can I distribute a Modelica package (under Modelica License 2) under a different license when I first add classes to the package?</h4>
<p>No. The package itself remains unmodified, i.e., it is Original Work, and therefore the license for this part must remain under Modelica License 2. The newly added classes can be, however, under a different license. </p>
<p><b>Can I copy a class out of a Modelica package (under Modelica License 2) and include it <u>unmodified</u> in a Modelica package under a <u>commercial/proprietary license</u>?</b></p>
<p>No, according to article 2c). However, you can include model, block, function, package, record and connector classes in your Modelica package under <u>Modelica License 2</u>. This means that your Modelica package could be under a commercial/proprietary license, but one or more classes of it are under Modelica License 2.</p>
<p>Note, a &ldquo;type&rdquo; class (e.g., type Angle = Real(unit=&rdquo;rad&rdquo;)) can be copied and included unmodified under a commercial/proprietary license (for details, see the next question). </p>
<p><b>Can I copy a type class or <u>part</u> of a model, block, function, record, connector class, out of a Modelica package (under Modelica License 2) and include it modified or unmodified in a Modelica package under a <u>commercial/proprietary</u> license</b></p>
<p>Yes, according to article 2d), since this will in the end usually qualify as Derivative Work. The reasoning is the following: A type class or part of another class (e.g., an equation, a declaration, part of a class description) cannot be utilized &ldquo;by its own&rdquo;. In order to make this &ldquo;usable&rdquo;, you have to add additional code in order that the class can be utilized. This is therefore usually Derivative Work and Derivative Work can be provided under a different license. Note, this only holds, if the additional code introduced is sufficient to qualify for Derivative Work. Merely, just copying a class and changing, say, one character in the documentation of this class would be no Derivative Work and therefore the copied code would have to stay under Modelica License 2. </p>
<p><b>Can I copy a class out of a Modelica package (under Modelica License 2) and include it in <u>modified </u>form in a <u>commercial/proprietary</u> Modelica package?</b></p>
<p>Yes. If the modification can be seen as a &ldquo;Derivative Work&rdquo;, you can place it under your commercial/proprietary license. If the modification does not qualify as &ldquo;Derivative Work&rdquo; (e.g., bug fixes, vendor specific annotations), it must remain under Modelica License 2. This means that your Modelica package could be under a commercial/proprietary license, but one or more parts of it are under Modelica License 2. </p>
<h4>Can I distribute a &ldquo;save total model&rdquo; under my commercial/proprietary license, even if classes under Modelica License 2 are included?</h4>
<p>Your classes of the &ldquo;save total model&rdquo; can be distributed under your commercial/proprietary license, but the classes under Modelica License 2 must remain under Modelica License 2. This means you can distribute a &ldquo;save total model&rdquo;, but some parts might be under Modelica License 2. </p>
<h4>Can I distribute a Modelica package (under Modelica License 2) in encrypted form?</h4>
<p>Yes. Note, if the encryption does not allow &ldquo;copying&rdquo; of classes (in to unencrypted Modelica source code), you have to send the Modelica source code of this package to your customer, if he/she wishes it, according to article&nbsp;6. </p>
<h4>Can I distribute an executable under my commercial/proprietary license, if the model from which the executable is generated uses models from a Modelica package under Modelica License 2?</h4>
<p>Yes, according to article 2d), since this is seen as Derivative Work. The reasoning is the following: An executable allows the simulation of a concrete model, whereas models from a Modelica package (without pre-processing, translation, tool run-time library) are not able to be simulated without tool support. By the processing of the tool and by its run-time libraries, significant new functionality is added (a model can be simulated whereas previously it could not be simulated) and functionality available in the package is removed (e.g., to build up a new model by dragging components of the package is no longer possible with the executable). </p>
<h4>Is my modification to a Modelica package (under Modelica License 2) a Derivative Work?</h4>
<p>It is not possible to give a general answer to it. To be regarded as &QUOT;an original work of authorship&QUOT;, a derivative work must be different enough from the original or must contain a substantial amount of new material. Making minor changes or additions of little substance to a preexisting work will not qualify the work as a new version for such purposes. </p>
<p><b></font><font style=\"color: #008000; \">Using or Distributing a Modelica <u>Document</u> under the Modelica License 2</font></b> </p>
<p><font style=\"font-size: 10pt; \">This section is devoted especially for the following applications:</p>
<p>A Modelica tool extracts information out of a Modelica package and presents the result in form of a &ldquo;manual&rdquo; for this package in, e.g., html, doc, or pdf format. </p>
<p>The Modelica language specification is a document defining the Modelica language. It will be licensed under Modelica License 2. </p>
<p>Someone writes a book about the Modelica language and/or Modelica packages and uses information which is available in the Modelica language specification and/or the corresponding Modelica package. </p>
<h4>Can I sell a manual that was basically derived by extracting information automatically from a Modelica package under Modelica License 2 (e.g., a &ldquo;reference guide&rdquo; of the Modelica Standard Library):</h4>
<p>Yes. Extracting information from a Modelica package, and providing it in a human readable, suitable format, like html, doc or pdf format, where the content is significantly modified (e.g. tables with interface information are constructed from the declarations of the public variables) qualifies as Derivative Work and there are no restrictions to charge a fee for Derivative Work under alternative 2d). </p>
<p><b>Can I copy a text passage out of a Modelica document (under Modelica License 2) and use it <u>unmodified</u> in my document (e.g. the Modelica syntax description in the Modelica Specification)?</b></p>
<p>Yes. In case you distribute your document, the copied parts are still under Modelica License 2 and you are not allowed to charge a license fee for this part. You can, of course, charge a fee for the rest of your document. </p>
<p><b>Can I copy a text passage out of a Modelica document (under Modelica License 2) and use it in <u>modified</u> form in my document?</b></p>
<p>Yes, the creation of Derivative Works is allowed. In case the content is significantly modified this qualifies as Derivative Work and there are no restrictions to charge a fee for Derivative Work under alternative 2d). </p>
<h4>Can I sell a printed version of a Modelica document (under Modelica License 2), e.g., the Modelica Language Specification?</h4>
<p>No, if you are not the copyright-holder, since article 2c) does not allow a selling fee for a (in this case physical) copy. However, mere printing and shipping costs may be recovered.</p>
</html>"));
    end ModelicaLicense2;

  package ReleaseNotes "Release notes"
    extends Modelica.Icons.ReleaseNotes;

  class Version_1_0 "Version 1.0.0 (28.4.2015)"
    extends Modelica.Icons.ReleaseNotes;

  annotation (Documentation(info="<html>
<ul>
<li>Separation the Chemical from Physiolibrary to https://github.com/MarekMatejak/Chemical from https://github.com/MarekMatejak/Physiolibrary branche PhysicalChemistry </li>
<li><font style=\"color: #333333; \">Components for solution, substance, chemical reaction, diffusion, gas dissolution, semipermeable membranes, chemical speciation of macromolecules, ..</font></li>
<li><font style=\"color: #333333; \">The library uses the Modelica Standard Libary (MSL) version 3.2.</font></li>
</ul>
</html>"));
  end Version_1_0;

  class Version_1_1 "Version 1.1.0 (20.5.2015)"
    extends Modelica.Icons.ReleaseNotes;

  annotation (Documentation(info="<html>
<ul>
<li>new state of matter - ideal gas</li>
<li>Solution with mechanical, thermal and electrical port </li>
<li>Ideal gas solution</li>
<li>Sensor of partial pressure</li>
<li>Sensor of dissociation constant of chemical reaction for hypothetical pure substances&rsquo; scheme </li>
<li>New Examples </li>
<li>New icon for electron transfer</li>
<li>New icon for chemical buffer </li>
</ul>
</html>"));
  end Version_1_1;
   annotation (Documentation(info="<html>
<p>This section summarizes the changes that have been performed on the Chemical. </p>
</html>"));

  end ReleaseNotes;

  class NewRelease "Publishing new release"
    extends Modelica.Icons.Information;

   annotation (Documentation(info="<html>
<p><br>New release must be numbered by Semantic Versioning 2.0.0, see <a href=\"http://semver.org/\">semver.org</a>. </p>
<p><br>If minor version, then the conversion script must be written and connected with package Chemical using &QUOT;annotation(conversion(from(version=..)))&QUOT;! </p>
<p><br>To clean the code from dummy annotations try to use script <a href=\"https://github.com/dietmarw/trimtrailingwhitespaces\">ttws</a>. </p>
<p>To check english spelling try to use <a href=\"https://github.com/vlajos/misspell_fixer\">missspell_fixer</a>.</p>
<p><br>Update version number to &QUOT;X.Y.Z&QUOT;: </p>
<ul>
<li>At package Chemical annotation: (version=&QUOT;X.Y.Z&QUOT;) together with &QUOT;versionBuild&QUOT;, &QUOT;versionDate&QUOT; and &QUOT;dateModified&QUOT; attribute </li>
<li>At file &QUOT;./Chemical/libraryinfo.mos&QUOT; </li>
</ul>
<p><br>Update release notes: </p>
<ul>
<li>At UsersGuide.ReleaseNotes</li>
<li>At file &QUOT;./README.md&QUOT;, together with update of &QUOT;Current release&QUOT; section.</li>
</ul>
<p><br>Publish release in GitHub: </p>
<ul>
<li>Prepare release in &QUOT;master&QUOT; branch</li>
<li>Install, Check, Test, Test, Test.. </li>
<li>Delete branch &QUOT;release&QUOT; </li>
<li>Create new branch &QUOT;release&QUOT; from &QUOT;master&QUOT; branch </li>
<li>Rename directory &QUOT;Chemical&QUOT; in release branch to directory &QUOT;Chemical X.Y.Z&QUOT; </li>
<li>Commint and Push </li>
<li>Draft a new release from &QUOT;release&QUOT; branch with number &QUOT;vX.Y.Z&QUOT; and with release notes. </li>
</ul>
</html>"));
  end NewRelease;

  class Contact "Contact"
    extends Modelica.Icons.Contact;

   annotation (Documentation(info="<html>
<p>Marek Matej&aacute;k</p>
<p>email: marek@matfy.cz</p>
<p>skype: marek.matejak</p>
<p>tel: +420 776 301 395</p>
</html>"));

  end Contact;

  annotation (__Dymola_DocumentationClass=true, Documentation(info="<html>
<p>Package <b>Chemical </b>is a modelica package for <b>Electro-Chemical processes </b>that is developed from <b>Physiolibrary</b> modelica implementation, see <a href=\"http://patf-biokyb.lf1.cuni.cz/wiki/hummod/hummod\">http://www.physiolibrary.org</a>. It provides connectors and model components fitted for electro-chemical models. </p>
</html>"));
  end UsersGuide;


 extends Modelica.Icons.Package;


  package Examples "Examples that demonstrate usage of chemical library"
  extends Modelica.Icons.ExamplesPackage;

    package Substances "Definitions of substances"
        extends Modelica.Icons.Package;

      constant Interfaces.Incompressible.SubstanceData Silver_solid(
        MolarWeight=0.1078682,
        z=0,
        DfH_25degC=0,
        DfG_25degC_1bar=0,
        Cp=25.4,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"}) "Ag(s)";

      constant Interfaces.Incompressible.SubstanceData Silver_aqueous(
        MolarWeight=0.1078682,
        z=1,
        DfH_25degC=105900,
        DfG_25degC_1bar=77100,
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Ag+(aq)";

      constant Interfaces.Incompressible.SubstanceData SilverChloride_solid(
        MolarWeight=0.14332,
        z=0,
        DfH_25degC=-127030,
        DfG_25degC_1bar=-109720,
        Cp=50.8,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "AgCl(s)";

      constant Interfaces.Incompressible.SubstanceData Calcium_aqueous(
        MolarWeight=0.0401,
        z=2,
        DfH_25degC=-542960,
        DfG_25degC_1bar=-542960 - 298.15*(33.67),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Ca++(aq)";

      constant Interfaces.Incompressible.SubstanceData Chloride_aqueous(
        MolarWeight=0.03545,
        z=-1,
        DfH_25degC=-167460,
        DfG_25degC_1bar=-131170,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "Cl-(aq)";

      constant Interfaces.IdealGas.SubstanceData CarbonDioxide_gas(
        MolarWeight=0.044,
        DfH_25degC=-393500,
        DfG_25degC_1bar=-394400,
        Cp=37.1,
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "CO2(g)";

      constant Interfaces.Incompressible.SubstanceData CarbonDioxide_aqueous(
        MolarWeight=0.044,
        DfH_25degC=-412900,
        DfG_25degC_1bar=-386200,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "CO2(aq)";

      constant Interfaces.Incompressible.SubstanceData Carbonate_aqueous(
        MolarWeight=0.06001,
        z=-2,
        DfH_25degC=-676300,
        DfG_25degC_1bar=-676300 - 298.15*(-497.065),
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "CO3--(aq)";

      constant Interfaces.Incompressible.SubstanceData Electrone_solid(
        MolarWeight=5.4857990946e-7,
        z=-1,
        DfH_25degC=0,
        DfG_25degC_1bar=0,
        Cp=0,
        References={"http://physics.nist.gov/cgi-bin/cuu/Value?mme, To solve standard electo-chemical cell potentials"}) "e-(s)";

      constant Interfaces.Incompressible.SubstanceData Iron2_aqueous(
        MolarWeight=0.05585,
        z=2,
        DfH_25degC=-87860,
        DfG_25degC_1bar=-87860 - 298.15*(-9.93),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Fe++(aq)";

      constant Interfaces.Incompressible.SubstanceData Iron3_aqueous(
        MolarWeight=0.05585,
        z=3,
        DfH_25degC=-47700,
        DfG_25degC_1bar=-47700 - 298.15*(-124.77),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Fe+++(aq)";

      constant Interfaces.Incompressible.SubstanceData Glucose_solid(
        MolarWeight=0.1806,
        DfH_25degC=-1274500,
        DfG_25degC_1bar=-1274500 - 298.15*(-1220.66),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "Glu(s)";

      constant Interfaces.IdealGas.SubstanceData Hydrogen_gas(
        MolarWeight=0.00201588,
        z=0,
        DfH_25degC=0,
        DfG_25degC_1bar=0,
        Cp=28.8,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"}) "H2(g)";

      constant Interfaces.Incompressible.SubstanceData CarbonicAcid_aqueous(
        MolarWeight=0.062027,
        DfH_25degC=-699700,
        DfG_25degC_1bar=-699700 - 298.15*(-256.582),
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H2CO3(aq)";

      constant Interfaces.IdealGas.SubstanceData Water_gas(
        MolarWeight=0.018015,
        DfH_25degC=-241830,
        DfG_25degC_1bar=-228590,
        Cp=33.6,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H2O(g)";

      constant Interfaces.Incompressible.SubstanceData Water_liquid(
        MolarWeight=0.018015,
        DfH_25degC=-285830,
        DfG_25degC_1bar=-237190,
        Cp=75.3,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H2O(l)";
     //   Cv=74.539,

      constant Interfaces.Incompressible.SubstanceData Water_IceIh(
        MolarWeight=0.018015,
        DfH_25degC=-292639,
        DfG_25degC_1bar=-236590,
        Cp=37.77,
        References={"http://www1.lsbu.ac.uk/water/water_properties.html#pot"})
      "H2O(s) - Ice I h";

      constant Interfaces.Incompressible.SubstanceData DihydrogenPhosphate_aqueous(
        MolarWeight=0.095,
        z=-1,
        DfH_25degC=-1302480,
        DfG_25degC_1bar=-1302480 - 298.15*(-561.395),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "H2PO4-(aq)";

      constant Interfaces.Incompressible.SubstanceData Hydronium_aqueous(
        MolarWeight=0.019022,
        z=1,
        DfH_25degC=-285840,
        DfG_25degC_1bar=-285840 - 298.15*(-163.17),
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H3O+(aq)";

      constant Interfaces.Incompressible.SubstanceData PhosphoricAcid_aqueous(
        MolarWeight=0.095,
        DfH_25degC=-1288000,
        DfG_25degC_1bar=-1288000 - 298.15*(-496.4),
        References={"https://en.wikipedia.org/wiki/Phosphoric_acid, https://www.researchgate.net/publication/6600409_Standard_thermodynamic_properties_of_H3PO4%28aq%29_over_a_wide_range_of_temperatures_and_pressures"})
      "H3PO4(aq)";

      constant Interfaces.Incompressible.SubstanceData Proton_aqueous(
        MolarWeight=0.001007,
        z=1,
        DfH_25degC=0,
        DfG_25degC_1bar=0,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "H+(aq)";
                 // as hypothetical HA <-> H+ + A- simplification of H2O + HA <-> H3O+ + A-";

      constant Interfaces.Incompressible.SubstanceData Bicarbonate_aqueous(
        MolarWeight=0.06102,
        z=-1,
        DfH_25degC=-691100,
        DfG_25degC_1bar=-691100 - 298.15*(-348.82),
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "HCO3-(aq)";

      constant Interfaces.Incompressible.SubstanceData Bicarbonate_blood(
        MolarWeight=0.06102,
        z=-1,
        DfH_25degC=-691100,
        DfG_25degC_1bar=-691100 - 298.15*(-348.82),
        gamma=0.79,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "HCO3-(blood)";

      constant Interfaces.Incompressible.SubstanceData HydrogenPhosphate_aqueous(
        MolarWeight=0.095,
        z=-2,
        DfH_25degC=-1298700,
        DfG_25degC_1bar=-1298700 - 298.15*(-686.232),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "HPO4--(aq)";

      constant Interfaces.Incompressible.SubstanceData HydrogenSulfate_aqueous(
        MolarWeight=0.097,
        z=-1,
        DfH_25degC=-885750,
        DfG_25degC_1bar=-752870,
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "HSO4-(aq)";

      constant Interfaces.Incompressible.SubstanceData Potassium_aqueous(
        MolarWeight=0.0391,
        z=1,
        DfH_25degC=-251200,
        DfG_25degC_1bar=-251200 - 298.15*(103.97),
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "K+(aq)";

      constant Interfaces.Incompressible.SubstanceData Magnesium_aqueous(
        MolarWeight=0.0243,
        z=2,
        DfH_25degC=-461960,
        DfG_25degC_1bar=-461960 - 298.15*(-19.99),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "Mg++(aq)";

      constant Interfaces.Incompressible.SubstanceData Sodium_aqueous(
        MolarWeight=0.02299,
        z=1,
        DfH_25degC=-239660,
        DfG_25degC_1bar=-239660 - 298.15*(74.49),
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "Na+(aq)";

      constant Interfaces.Incompressible.SubstanceData Amonium_aqueous(
        MolarWeight=0.01804,
        z=1,
        DfH_25degC=-132800,
        DfG_25degC_1bar=-132800 - 298.15*(-178.77),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "NH4+(aq)";

      constant Interfaces.IdealGas.SubstanceData Oxygen_gas(
        MolarWeight=0.032,
        DfH_25degC=0,
        DfG_25degC_1bar=0,
        Cp=29.4,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"}) "O2(g)";

      constant Interfaces.Incompressible.SubstanceData Oxygen_aqueous(
        MolarWeight=0.032,
        DfH_25degC=-11700,
        DfG_25degC_1bar=16320,
        References={"http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.pdf, https://books.google.cz/books?id=dr-VBAAAQBAJ&pg=PA156&lpg=PA156&dq=Gibbs+energy+of+formation++%22O2(aq)%22&source=bl&ots=09N5CxY7OD&sig=hbsTXQvX59vXBqHUjFVVIZQpHCA&hl=cs&sa=X&ei=sDQtVaeUMMaRsAHpzYHgAg&redir_esc=y#v=onepage&q=Gibbs%20energy%20of%20formation%20%20%22O2(aq)%22&f=false"})
      "O2(aq)";

      constant Interfaces.Incompressible.SubstanceData Hydroxide_aqueous(
        MolarWeight=0.017006,
        z=-1,
        DfH_25degC=-229940,
        DfG_25degC_1bar=-157300,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html"})
      "OH-(aq)";

      constant Interfaces.Incompressible.SubstanceData Lead_solid(
        MolarWeight=0.2072,
        z=0,
        DfH_25degC=0,
        DfG_25degC_1bar=0,
        Cp=26.4,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"}) "Pb(s)";

      constant Interfaces.Incompressible.SubstanceData LeadDioxide_solid(
        MolarWeight=0.2391988,
        z=0,
        DfH_25degC=-276600,
        DfG_25degC_1bar=-219000,
        Cp=64.6,
        References={"http://www.vias.org/genchem/standard_enthalpies_table.html, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"})
      "PbO2(s)";

      constant Interfaces.Incompressible.SubstanceData LeadSulfate_solid(
        MolarWeight=0.30326,
        z=0,
        DfH_25degC=-918400,
        DfG_25degC_1bar=-811200,
        Cp=103.2,
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf"})
      "PbSO4(s)";

      constant Interfaces.Incompressible.SubstanceData Phosphate_aqueous(
        MolarWeight=0.095,
        z=-3,
        DfH_25degC=-1284070,
        DfG_25degC_1bar=-1284070 - 298.15*(-866.946),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "PO4---(aq)";

      constant Interfaces.Incompressible.SubstanceData Sulphates_aqueous(
        MolarWeight=0.09607,
        z=-2,
        DfH_25degC=-907500,
        DfG_25degC_1bar=-907500 - 298.15*(-555.123),
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf"})
      "SO4--(aq)";

      constant Interfaces.Incompressible.SubstanceData Ethanol_liquid(
        MolarWeight=0.04607,
        z=0,
        DfH_25degC=-276980,
        DfG_25degC_1bar=-174180,
        Cp=112.4,
        density=789,
        References={"http://www.mhhe.com/physsci/chemistry/chang7/ssg/graphics/chang7/pdf/cng7pa08.pdf, https://en.wikipedia.org/wiki/Ethanol_(data_page)"})
      "Ethanol C2H5OH(l)";

        //Some organic molecules: https://www.e-education.psu.edu/drupal6/files/be497b/pdf/Bioenergetics_AppA.pdf
        //Other source: http://www.update.uu.se/~jolkkonen/pdf/CRC_TD.pdf
    end Substances;

    model SimpleReaction
    "The simple chemical reaction A<->B with equilibrium B/A = 2"
       extends Modelica.Icons.Example;

      constant Real K = 2 "Dissociation constant of the reaction";

      constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

    Components.SimpleSolution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Components.Substance A(
        amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-52,-8},{-32,12}})));

      Components.Reaction reaction annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
      Components.Substance B(
        substanceData( DfG_25degC_1bar=-R*T_25degC*log(K)),
        amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{62,-8},{42,12}})));

    equation
      connect(reaction.products[1], B.port_a) annotation (Line(
          points={{10,2},{42,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(A.solution, solution.solution) annotation (Line(
          points={{-48,-8},{-48,-92},{0,-92},{0,-98}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(B.solution, solution.solution) annotation (Line(points={{58,-8},{
            58,-92},{0,-92},{0,-98}},    smooth=Smooth.None));
      connect(A.port_a, reaction.substrates[1]) annotation (Line(
          points={{-32,2},{-10,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
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

      constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

    Components.SimpleSolution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Components.Substance A(amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
      Components.Reaction reaction(nS=2)
        annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Components.Substance B(amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Components.Substance C(substanceData(DfG_25degC_1bar=-R*T_25degC*log(Kx)), amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{68,-8},{48,12}})));

    equation
      connect(reaction.products[1], C.port_a) annotation (Line(
          points={{24,2},{48,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(A.solution, solution.solution) annotation (Line(
          points={{-30,2},{-30,-90},{0,-90},{0,-98}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(C.solution, solution.solution) annotation (Line(points={{64,-8},{
            66,-8},{66,-90},{0,-90},{0,-98}},    smooth=Smooth.None));
      connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},
            {-30,-90},{0,-90},{0,-98}},  smooth=Smooth.None));

      connect(B.port_a, reaction.substrates[1]) annotation (Line(
          points={{-14,-14},{-10,-14},{-10,0},{4,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(A.port_a, reaction.substrates[2]) annotation (Line(
          points={{-14,12},{-10,12},{-10,4},{4,4}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.001));
    end SimpleReaction2;

    model HeatingOfWater "Heating of 1 kg water"
      extends Modelica.Icons.Example;

    Components.Solution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
      Components.Substance H2O(
        redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible,
        substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=55.508)
        annotation (Placement(transformation(extent={{-4,-16},{16,4}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
        annotation (Placement(transformation(extent={{28,-76},{48,-56}})));
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{-70,60},{-50,80}})));
    equation
      connect(H2O.solution, solution.solution) annotation (Line(
          points={{0,-16},{0,-98}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(fixedHeatFlow.port, solution.heatPort) annotation (Line(
          points={{48,-66},{60,-66},{60,-98}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(ground.p, solution.electricPin) annotation (Line(
          points={{-60,80},{-60,100}},
          color={0,0,255},
          smooth=Smooth.None));
      annotation (experiment(StopTime=1),
      Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end HeatingOfWater;

    model HeatingOfAlcohol "Heating of 50% ethanol"
      extends Modelica.Icons.Example;

    Components.Solution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
                                /*(mass_start=0.5 + (55.508/2)*Substances.Ethanol_liquid.MolarWeight,
    volume_start=1/(0.997*0.91251))*/
      Components.Substance H2O(
        redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible,
        substanceData=Chemical.Examples.Substances.Water_liquid,
      amountOfSubstance_start=55.508/2)
        annotation (Placement(transformation(extent={{-46,-8},{-26,12}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
        annotation (Placement(transformation(extent={{28,-76},{48,-56}})));
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{-70,62},{-50,82}})));
      Components.Substance Ethanol(
      redeclare package stateOfMatter =
          Chemical.Interfaces.Incompressible,
      substanceData=Chemical.Examples.Substances.Ethanol_liquid,
      amountOfSubstance_start=55.508/2)
      annotation (Placement(transformation(extent={{18,-8},{38,12}})));
    equation
      connect(H2O.solution, solution.solution) annotation (Line(
          points={{-42,-8},{-42,-34},{0,-34},{0,-98}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(fixedHeatFlow.port, solution.heatPort) annotation (Line(
          points={{48,-66},{60,-66},{60,-98}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(ground.p, solution.electricPin) annotation (Line(
          points={{-60,82},{-60,100}},
          color={0,0,255},
          smooth=Smooth.None));
    connect(solution.solution, Ethanol.solution) annotation (Line(
        points={{0,-98},{0,-34},{22,-34},{22,-8}},
        color={158,66,200},
        smooth=Smooth.None));
      annotation (experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end HeatingOfAlcohol;

    model ExothermicReaction
    "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"

       extends Modelica.Icons.Example;

      parameter Modelica.SIunits.MolarEnergy ReactionEnthalpy=-55000;

    Components.Solution thermal_isolated_solution
      annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));
      Components.Substance A( amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
      Components.Reaction reaction
        annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
      Components.Substance B( amountOfSubstance_start=0.1, substanceData(DfH_25degC=ReactionEnthalpy))
        annotation (Placement(transformation(extent={{40,-60},{20,-40}})));

    Components.Solution solution_at_constant_temperature
      annotation (Placement(transformation(extent={{-100,0},{98,94}})));
      Components.Substance A1(amountOfSubstance_start=0.9)
        annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      Components.Reaction reaction1
        annotation (Placement(transformation(extent={{-8,40},{12,60}})));
      Components.Substance B1(amountOfSubstance_start=0.1, substanceData(DfH_25degC=ReactionEnthalpy))
        annotation (Placement(transformation(extent={{40,40},{20,60}})));

      Modelica.SIunits.HeatFlowRate q
      "Heat flow to environment to reach constant temperature";
      Modelica.SIunits.Temperature t
      "Temperature if the solution is ideally thermal isolated from environment";
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{-58,-42},{-38,-22}})));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=
          298.15) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          origin={82,50})));
      Components.Substance H2O(
        redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible,
        substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=55.508)
        annotation (Placement(transformation(extent={{20,4},{40,24}})));
      Components.Substance H2O1(
        redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible,
        substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=55.508)
        annotation (Placement(transformation(extent={{20,-94},{40,-74}})));
    equation
      q = -solution_at_constant_temperature.heatPort.Q_flow;

      t = thermal_isolated_solution.solution.T;

      connect(A.port_a, reaction.substrates[1]) annotation (Line(
          points={{-20,-50},{-8,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(reaction.products[1], B.port_a) annotation (Line(
          points={{12,-50},{20,-50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(B.solution, thermal_isolated_solution.solution) annotation (Line(
          points={{36,-60},{36,-64},{-1,-64},{-1,-99.06}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(A.solution, thermal_isolated_solution.solution) annotation (Line(
            points={{-36,-60},{-36,-64},{-1,-64},{-1,-99.06}},
                                                             smooth=Smooth.None));
      connect(A1.port_a, reaction1.substrates[1]) annotation (Line(
          points={{-20,50},{-8,50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(reaction1.products[1], B1.port_a) annotation (Line(
          points={{12,50},{20,50}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(B1.solution, solution_at_constant_temperature.solution) annotation (
          Line(
          points={{36,40},{36,34},{-1,34},{-1,0.94}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(A1.solution, solution_at_constant_temperature.solution) annotation (
          Line(points={{-36,40},{-36,34},{-1,34},{-1,0.94}},
                                                          smooth=Smooth.None));
    connect(solution_at_constant_temperature.electricPin, ground.p) annotation (
       Line(
        points={{-60.4,94},{-60,94},{-60,-22},{-48,-22}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(thermal_isolated_solution.electricPin, ground.p) annotation (Line(
        points={{-60.4,-6},{-60,-6},{-60,-22},{-48,-22}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(fixedTemperature.port, solution_at_constant_temperature.heatPort)
      annotation (Line(
        points={{72,50},{58.4,50},{58.4,0.94}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(solution_at_constant_temperature.solution, H2O.solution)
      annotation (Line(
        points={{-1,0.94},{24,0.94},{24,4}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(thermal_isolated_solution.solution, H2O1.solution) annotation (Line(
        points={{-1,-99.06},{24,-99.06},{24,-94}},
        color={158,66,200},
        smooth=Smooth.None));
      annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=0.001));
    end ExothermicReaction;

    model PowerGeneration "Hydrogen burning piston"
      extends Modelica.Icons.Example;

      parameter Modelica.SIunits.Volume V=0.001 "Initial volume";
     // parameter Modelica.SIunits.Pressure p=100000 "Initial pressure";
      parameter Modelica.SIunits.Temperature T=298.15 "Initial temperature";

      parameter Modelica.SIunits.Area A=0.01 "Cross area of cylinder";

      //p*V=n*R*T
     // parameter Modelica.SIunits.AmountOfSubstance n=p*V/(Modelica.Constants.R*T)
     //   "Initial amount of substances in sulution";

      Components.IdealGasSolution idealGas(
        SurfaceArea=A)
        annotation (Placement(transformation(extent={{-50,-68},{50,32}})));
                       // AmbientPressure=p)
      //  volume_start=V,
      Components.Substance H2_gas(
      redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Examples.Substances.Hydrogen_gas,
        amountOfSubstance_start(displayUnit="mmol") = 0.026)
      annotation (Placement(transformation(extent={{-40,-44},{-20,-24}})));
      Components.Substance O2_gas(
      substanceData=Chemical.Examples.Substances.Oxygen_gas,
      redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGas,
        amountOfSubstance_start(displayUnit="mmol") = 0.013)
      annotation (Placement(transformation(extent={{-40,-8},{-20,12}})));
      Components.Substance H2O_gas(substanceData=Chemical.Examples.Substances.Water_gas,
        redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGas)
      annotation (Placement(transformation(extent={{44,-26},{24,-6}})));
      Components.Reaction reaction(
      nS=2,
      s={2,1},
      p={2}) annotation (Placement(transformation(extent={{-10,-26},{10,-6}})));
      Modelica.Mechanics.Translational.Components.Spring spring(c=1e6) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={0,52})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=2)
        annotation (Placement(transformation(extent={{44,-96},{64,-76}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature coolerTemperature(T=298.15)
        annotation (Placement(transformation(extent={{96,-96},{76,-76}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed
        annotation (Placement(transformation(extent={{8,72},{28,92}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed1
        annotation (Placement(transformation(extent={{-10,-98},{10,-78}})));
    equation
    connect(reaction.products[1], H2O_gas.port_a) annotation (Line(
        points={{10,-16},{24,-16}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(H2_gas.port_a, reaction.substrates[1]) annotation (Line(
        points={{-20,-34},{-16,-34},{-16,-18},{-10,-18}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(O2_gas.port_a, reaction.substrates[2]) annotation (Line(
        points={{-20,2},{-16,2},{-16,-14},{-10,-14}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(H2_gas.solution, idealGas.solution) annotation (Line(
        points={{-36,-44},{-44,-44},{-44,-67},{30,-67}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(O2_gas.solution, idealGas.solution) annotation (Line(
        points={{-36,-8},{-44,-8},{-44,-67},{30,-67}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(H2O_gas.solution, idealGas.solution) annotation (Line(
        points={{40,-26},{40,-48},{30,-48},{30,-67}},
        color={158,66,200},
        smooth=Smooth.None));
      connect(idealGas.surfaceFlange, spring.flange_a) annotation (Line(
          points={{0,32},{0,42}},
          color={0,127,0},
          smooth=Smooth.None));
      connect(idealGas.heatPort, thermalConductor.port_a) annotation (Line(
          points={{-30,-69},{-30,-86},{44,-86}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(thermalConductor.port_b, coolerTemperature.port) annotation (Line(
          points={{64,-86},{76,-86}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(fixed.flange, spring.flange_b) annotation (Line(
          points={{18,82},{18,88},{0,88},{0,62},{6.66134e-016,62}},
          color={0,127,0},
          smooth=Smooth.None));
    connect(idealGas.bottom, fixed1.flange) annotation (Line(
        points={{0,-69},{0,-88}},
        color={0,127,0},
        smooth=Smooth.None));
      annotation ( experiment(StopTime=1), Diagram(graphics));
    end PowerGeneration;

    model WaterVaporization "Evaporation of water"
       extends Modelica.Icons.Example;

       parameter Modelica.SIunits.Temperature T_start=273.15
      "Initial temperature";

    Components.Solution liquid(temperature_start=T_start)
        annotation (Placement(transformation(extent={{-98,-98},{-6,-8}})));

      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,

      Components.Substance H2O_liquid(substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=55.508) "Liquid water"
        annotation (Placement(transformation(extent={{-54,-64},{-74,-44}})));

    Components.IdealGasSolution gas(temperature_start=T_start)
        annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                    /*volume_start(
        displayUnit="l") = 0.001, */
      Components.GasSolubility gasSolubility(useWaterCorrection=false, KC=10)
        annotation (Placement(transformation(extent={{-98,24},{-78,44}})));
      Components.Substance H2O_gaseuous(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas,               substanceData=Chemical.Examples.Substances.Water_gas,
      amountOfSubstance_start(displayUnit="mmol") = 0.001)
        annotation (Placement(transformation(extent={{-4,54},{-24,74}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                           fixedTemperature
                  annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          origin={84,8})));
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{-80,-36},{-60,-16}})));
      Modelica.Blocks.Sources.Clock clock(offset=1*T_start)
        annotation (Placement(transformation(extent={{62,36},{82,56}})));
    Components.Substance otherSubstances(substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-36,14},{-16,34}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(
        G=1e6) annotation (Placement(transformation(extent={{48,-6},{68,14}})));
    equation

      connect(H2O_liquid.solution, liquid.solution) annotation (Line(
          points={{-58,-64},{-58,-97.1},{-52,-97.1}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(H2O_liquid.port_a, gasSolubility.liquid_port) annotation (Line(
          points={{-74,-54},{-88,-54},{-88,24}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(gas.solution, H2O_gaseuous.solution) annotation (Line(
          points={{0,6.9},{-8,6.9},{-8,54}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H2O_gaseuous.port_a, gasSolubility.gas_port) annotation (Line(
          points={{-24,64},{-88,64},{-88,44}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(liquid.electricPin, ground.p) annotation (Line(
          points={{-79.6,-8},{-79.6,-8},{-70,-8},{-70,-16}},
          color={0,0,255},
          smooth=Smooth.None));
    connect(fixedTemperature.T, clock.y) annotation (Line(
        points={{96,8},{98,8},{98,46},{83,46}},
        color={0,0,127},
        smooth=Smooth.Bezier));
    connect(gas.solution, otherSubstances.solution) annotation (Line(
        points={{0,6.9},{-24,6.9},{-24,14},{-32,14}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(fixedTemperature.port, thermalConductor.port_b) annotation (Line(
        points={{74,8},{72,8},{72,4},{68,4}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(gas.heatPort, thermalConductor.port_a) annotation (Line(
        points={{27.6,6.9},{34,6.9},{34,4},{48,4}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(liquid.heatPort, thermalConductor.port_a) annotation (Line(
        points={{-24.4,-97.1},{-2,-97.1},{-2,4},{48,4}},
        color={191,0,0},
        smooth=Smooth.None));
      annotation (
        experiment(StopTime=100),
        Documentation(info="<html>
<p>Please note, that the total content of CO2 and O2 in blood plasma and erythrocytes must be determined by including bicarbonate and hemoglobin connected amounts.  </p>
</html>", revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end WaterVaporization;

    model WaterSublimation "Sublimation of water"
       extends Modelica.Icons.Example;

       parameter Modelica.SIunits.Temperature T_start=273.15-50
      "Initial temperature";

      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,

    Components.IdealGasSolution gas(temperature_start=T_start, AmbientPressure=
          600)
        annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                    /*volume_start(
        displayUnit="l") = 0.001, */
      Components.Substance H2O_gaseuous(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas,               substanceData=Chemical.Examples.Substances.Water_gas,
      amountOfSubstance_start(displayUnit="mmol") = 0.001)
        annotation (Placement(transformation(extent={{-4,54},{-24,74}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                           fixedTemperature
                  annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          origin={84,8})));
      Modelica.Blocks.Sources.Clock clock(offset=1*T_start)
        annotation (Placement(transformation(extent={{62,36},{82,56}})));
    Components.Substance otherSubstances(substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=1)
      annotation (Placement(transformation(extent={{-36,14},{-16,34}})));
    Components.Solution solid(temperature_start=T_start, AmbientPressure=600)
        annotation (Placement(transformation(extent={{8,-98},{100,-8}})));
      Components.Substance H2O_solid(amountOfSubstance_start=55.508, substanceData=
            Chemical.Examples.Substances.Water_IceIh) "Solid water"
        annotation (Placement(transformation(extent={{70,-62},{50,-42}})));
      Components.GasSolubility gasSolubility1(
                                             useWaterCorrection=false, KC=10)
        annotation (Placement(transformation(extent={{-76,18},{-56,38}})));
      Modelica.Electrical.Analog.Basic.Ground ground1
        annotation (Placement(transformation(extent={{12,-38},{32,-18}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=1e6)
        annotation (Placement(transformation(extent={{48,-6},{68,14}})));
    equation

      connect(gas.solution, H2O_gaseuous.solution) annotation (Line(
          points={{0,6.9},{-8,6.9},{-8,54}},
          color={158,66,200},
          smooth=Smooth.None));
    connect(fixedTemperature.T, clock.y) annotation (Line(
        points={{96,8},{98,8},{98,46},{83,46}},
        color={0,0,127},
        smooth=Smooth.Bezier));
    connect(gas.solution, otherSubstances.solution) annotation (Line(
        points={{0,6.9},{-24,6.9},{-24,14},{-32,14}},
        color={158,66,200},
        smooth=Smooth.None));
      connect(solid.solution, H2O_solid.solution) annotation (Line(
          points={{54,-97.1},{60,-97.1},{60,-62},{66,-62}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H2O_gaseuous.port_a, gasSolubility1.gas_port) annotation (Line(
          points={{-24,64},{-66,64},{-66,38}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(gasSolubility1.liquid_port, H2O_solid.port_a) annotation (Line(
          points={{-66,18},{-66,-52},{50,-52}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(solid.electricPin, ground1.p) annotation (Line(
          points={{26.4,-8},{22,-8},{22,-18}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(fixedTemperature.port, thermalConductor.port_b) annotation (Line(
          points={{74,8},{72,8},{72,4},{68,4}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(gas.heatPort, thermalConductor.port_a) annotation (Line(
          points={{27.6,6.9},{34,6.9},{34,4},{48,4}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(solid.heatPort, thermalConductor.port_a) annotation (Line(
          points={{81.6,-97.1},{-2,-97.1},{-2,4},{48,4}},
          color={191,0,0},
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=50.01),
        Documentation(info="<html>
<p>Please note, that the total content of CO2 and O2 in blood plasma and erythrocytes must be determined by including bicarbonate and hemoglobin connected amounts.  </p>
</html>", revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end WaterSublimation;

    model GasSolubility "Dissolution of gases in liquids"
       extends Modelica.Icons.Example;

    Components.SimpleSolution blood_plasma
      annotation (Placement(transformation(extent={{-100,-76},{-8,14}})));
                                          //(amountOfSolution_start=52.3)
    Components.SimpleSolution red_cells
      annotation (Placement(transformation(extent={{8,-78},{102,14}})));
                                       //(amountOfSolution_start=39.7)

      Components.GasSolubility CO2_dissolutionP
      annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,

      Components.Substance CO2_unbound_plasma(substanceData=Substances.CarbonDioxide_aqueous)
      "Free dissolved CO2 in blood plasma"
      annotation (Placement(transformation(extent={{-90,-24},{-70,-4}})));
      Components.GasSolubility O2_dissolutionP
      annotation (Placement(transformation(extent={{-34,44},{-14,64}})));

    Sources.ExternalIdealGasSubstance O2_g_n1(
      substanceData=Substances.Oxygen_gas,
      PartialPressure=12665.626804425,
      TotalPressure=101325.0144354)
      annotation (Placement(transformation(extent={{22,76},{42,96}})));
      Components.Substance O2_unbound_plasma(substanceData=Substances.Oxygen_aqueous)
      "Free dissolved O2 in blood plasma"
      annotation (Placement(transformation(extent={{-50,-26},{-30,-6}})));
      Components.GasSolubility CO2_dissolutionE
      annotation (Placement(transformation(extent={{36,44},{56,64}})));

    Sources.ExternalIdealGasSubstance CO2_g_n2(
      substanceData=Substances.CarbonDioxide_gas,
      PartialPressure=5332.8954966,
      TotalPressure=101325.0144354)
      annotation (Placement(transformation(extent={{-56,78},{-36,98}})));

      Components.Substance CO2_unbound_erythrocyte(substanceData=Substances.CarbonDioxide_aqueous)
      "Free dissolved CO2 in red cells"
      annotation (Placement(transformation(extent={{18,-32},{38,-12}})));

      Components.GasSolubility O2_dissolutionE_NIST(useWaterCorrection=true)
      annotation (Placement(transformation(extent={{78,44},{98,64}})));
      Components.Substance O2_unbound_erythrocyte_NIST(substanceData=Chemical.Examples.Substances.Oxygen_aqueous)
      "Free dissolved O2 in red cells"
      annotation (Placement(transformation(extent={{58,-32},{78,-12}})));
    Components.Substance otherSubstances(substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=52.3)
      annotation (Placement(transformation(extent={{-42,-70},{-22,-50}})));
    Components.Substance otherSubstances_erythrocytes(substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=38.7)
      annotation (Placement(transformation(extent={{64,-68},{84,-48}})));
    equation

    connect(CO2_g_n2.port_a, CO2_dissolutionP.gas_port) annotation (Line(
        points={{-36,88},{-26,88},{-26,72},{-68,72},{-68,64}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(CO2_g_n2.port_a, CO2_dissolutionE.gas_port) annotation (Line(
        points={{-36,88},{-26,88},{-26,72},{46,72},{46,64}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(CO2_dissolutionP.liquid_port, CO2_unbound_plasma.port_a)
      annotation (Line(
        points={{-68,44},{-68,-14},{-70,-14}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(CO2_dissolutionE.liquid_port, CO2_unbound_erythrocyte.port_a)
      annotation (Line(
        points={{46,44},{46,-22},{38,-22}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(O2_g_n1.port_a, O2_dissolutionP.gas_port) annotation (Line(
        points={{42,86},{66,86},{66,70},{-24,70},{-24,64}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(O2_dissolutionP.liquid_port, O2_unbound_plasma.port_a) annotation (
        Line(
        points={{-24,44},{-24,-16},{-30,-16}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(CO2_unbound_plasma.solution, blood_plasma.solution) annotation (
        Line(
        points={{-86,-24},{-86,-75.1},{-54,-75.1}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(O2_unbound_plasma.solution, blood_plasma.solution) annotation (Line(
        points={{-46,-26},{-46,-75.1},{-54,-75.1}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(CO2_unbound_erythrocyte.solution, red_cells.solution) annotation (
        Line(
        points={{22,-32},{22,-77.08},{55,-77.08}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(O2_g_n1.port_a, O2_dissolutionE_NIST.gas_port) annotation (Line(
        points={{42,86},{66,86},{66,70},{88,70},{88,64}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(O2_dissolutionE_NIST.liquid_port, O2_unbound_erythrocyte_NIST.port_a)
      annotation (Line(
        points={{88,44},{88,-22},{78,-22}},
        color={158,66,200},
        thickness=1,
        smooth=Smooth.None));
    connect(O2_unbound_erythrocyte_NIST.solution, red_cells.solution)
      annotation (Line(
        points={{62,-32},{62,-77.08},{55,-77.08}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(blood_plasma.solution, otherSubstances.solution) annotation (Line(
        points={{-54,-75.1},{-38,-75.1},{-38,-70}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(red_cells.solution, otherSubstances_erythrocytes.solution)
      annotation (Line(
        points={{55,-77.08},{68,-77.08},{68,-68}},
        color={158,66,200},
        smooth=Smooth.None));
      annotation (
        experiment(StopTime=1),
        Documentation(info="<html>
<p>Please note, that the total content of CO2 and O2 in blood plasma and erythrocytes must be determined by including bicarbonate and hemoglobin connected amounts.  </p>
</html>", revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end GasSolubility;

    model EnzymeKinetics "Basic enzyme kinetics"
      extends Modelica.Icons.Example;

    Components.SimpleSolution solution
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
      Components.Substance    P(             substanceData(DfG_25degC_1bar=-Modelica.Constants.R*298.15*50))
        annotation (Placement(transformation(extent={{92,-12},{72,8}})));
      Components.Substance    S(amountOfSubstance_start=1)
        annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));

         parameter Modelica.SIunits.AmountOfSubstance tE=1e-6
      "Total amount of enzyme";
         parameter Real k_cat(unit="1/s", displayUnit="1/min")= 1
      "Forward rate of second reaction";
         constant Modelica.SIunits.Concentration Km=0.1
      "Michaelis constant = substrate concentration at rate of half Vmax";

        parameter Modelica.SIunits.MolarFlowRate Vmax=1e-5*k_cat
      "Maximal molar flow";
        parameter Modelica.SIunits.AmountOfSubstance AmountOfSolution= 55.508
      "Amount of solution used in kinetics";

          Components.Substance ES(substanceData(DfG_25degC_1bar=-Modelica.Constants.R*298.15*log(2/Km)),
          amountOfSubstance_start=tE/2)
            annotation (Placement(transformation(extent={{-12,-10},{8,10}})));
          Components.Substance E(amountOfSubstance_start=tE/2)
            annotation (Placement(transformation(extent={{-10,38},{10,58}})));
      Components.Reaction chemicalReaction(nS=2, KC=Vmax/(2*Modelica.Constants.R*298.15*log(2)))
        annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));

      Components.Reaction chemicalReaction1(nP=2, KC=Vmax/(2*Modelica.Constants.R*298.15*(50 - log(2))))
        annotation (Placement(transformation(extent={{24,-10},{44,10}})));

    Components.Substance otherSubstances(substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=52.3)
      annotation (Placement(transformation(extent={{42,-76},{62,-56}})));
    equation
         //Michaelis-Menton: v=((E.q_out.conc + ES.q_out.conc)*k_cat)*S.concentration/(Km+S.concentration);

      connect(S.port_a, chemicalReaction.substrates[1]) annotation (Line(
          points={{-72,-4},{-56,-4},{-56,-2},{-42,-2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(chemicalReaction.products[1], ES.port_a) annotation (Line(
          points={{-22,0},{8,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(ES.port_a, chemicalReaction1.substrates[1]) annotation (Line(
          points={{8,0},{24,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(E.port_a, chemicalReaction.substrates[2]) annotation (Line(
          points={{10,48},{-52,48},{-52,2},{-42,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(E.port_a, chemicalReaction1.products[2]) annotation (Line(
          points={{10,48},{54,48},{54,2},{44,2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(chemicalReaction1.products[1], P.port_a) annotation (Line(
          points={{44,-2},{58,-2},{58,-2},{72,-2}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(E.solution, solution.solution) annotation (Line(
          points={{-6,38},{-8,38},{-8,-98},{0,-98}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(ES.solution, solution.solution)
        annotation (Line(points={{-8,-10},{-8,-98},{0,-98}},          smooth=Smooth.None));
      connect(S.solution, solution.solution) annotation (Line(
          points={{-88,-14},{-88,-98},{0,-98}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(P.solution, solution.solution) annotation (Line(
          points={{88,-12},{88,-98},{0,-98}},
          color={158,66,200},
          smooth=Smooth.None));
    connect(solution.solution, otherSubstances.solution) annotation (Line(
        points={{0,-98},{46,-98},{46,-76}},
        color={158,66,200},
        smooth=Smooth.None));
          annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Be carefull, the assumption for Michaelis-Menton are very strong: </p>
<p>The substrate must be in sufficiently high concentration and the product must be in very low concentration to reach almost all enzyme in enzyme-substrate complex all time. ([S] &GT;&GT; Km) &AMP;&AMP; ([P] &LT;&LT; K2)</p>
<p><br>To recalculate the enzyme kinetics from Michaelis-Menton parameters Km, tE a k_cat is selected the same half-rate of the reaction defined as:</p>
<p>E = ES = tE/2 .. the amount of free enzyme is the same as the amount of enzyme-substrate complexes</p>
<p>S = Km .. the amount of substrate is Km</p>
<p>r = Vmax/2 = tE*k_cat / 2 .. the rate of reaction is the half of maximal rate</p>
<p><br>Conversions of molar concentration to mole fraction (MM is molar mass of the solvent in solution -&GT; 55.508 kg/mol for water):</p>
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
<p>-kC2 * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) ?=&LT;? Vmax * (1 + 1/( -uP&deg;/(R*T*ln(2)) - 1) )</p>
<p><br>(Vmax/2) * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) / ( - uP&deg; - R * T * ln(2) ) ?=&LT;? Vmax*(1 + R*T*ln(2) / ( -uP&deg; - R*T*ln(2)) )</p>
<p>(uP&deg; +<b> </b>R*T*ln(2/x(Km)) + R*T*ln(xP*xE/xES) ) ?=&LT;? 2*( - uP&deg; - R * T * ln(2) ) + 2*R*T*ln(2)</p>
<p>R*T*ln(xP*xE/xES) ?=&LT;? - uP&deg; - R*T*ln(2/x(Km)) </p>
<p>xP*xE/xES ?=&LT;? exp((- uP&deg; - R*T*ln(2/x(Km))/(R*T))</p>
<p>The equality is the equation of the equilibrium: xP*xE/xES = exp((- uP&deg; - uE&deg; + uES&deg; )/(R*T)) = exp((- uP&deg; - R*T*ln(2/x(Km))/(R*T))</p>
<p>If the equilibrium of the reaction is reached only by forward rate then xP*xE/xES must be less than the dissociation constant.</p>
<h4>The increasing of the amount of the enzyme</h4>
<p>In the situation of doubled amount of enzyme should double also the maximal speed of the reaction, shouldn&apos;t?</p>
<p>The assumptions of</p>
</html>"),
        experiment(StopTime=200000));
    end EnzymeKinetics;

    model ElectrochemicalCell
    "The electrochemical cell: Pt(s) | H2(g) | H+(aq), Cl-(aq) | AgCl(s) | Ag(s)"
     extends Modelica.Icons.Example;

    Components.SimpleSolution cathode(ElectricGround=false)
      annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
    Components.SimpleSolution anode(ElectricGround=false)
      annotation (Placement(transformation(extent={{62,-50},{96,50}})));

    Components.SimpleSolution solution1(ElectricGround=false)
      annotation (Placement(transformation(extent={{-30,-60},{38,6}})));

      Components.Substance  Ag(substanceData=
            Chemical.Examples.Substances.Silver_solid, amountOfSubstance_start=
          1)
        annotation (Placement(transformation(extent={{-72,-30},{-52,-10}})));
      Components.Substance Cl(substanceData=
            Chemical.Examples.Substances.Chloride_aqueous,
        amountOfSubstance_start=1)     annotation (Placement(transformation(extent={{-2,-36},
              {-22,-16}})));
      Components.Substance  AgCl(substanceData=
            Chemical.Examples.Substances.SilverChloride_solid)
        annotation (Placement(transformation(extent={{-76,4},{-56,24}})));
    Sources.ExternalIdealGasSubstance H2(
      substanceData=Chemical.Examples.Substances.Hydrogen_gas,
      PartialPressure=100000,
      TotalPressure=100000)
      annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Components.Substance H(substanceData=
            Chemical.Examples.Substances.Proton_aqueous,
        amountOfSubstance_start=1)    annotation (Placement(transformation(extent={{6,-36},
              {26,-16}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Components.Reaction electrodeReaction(nP=2, p={2,2}) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={52,6})));
      Components.Reaction electrodeReaction1(nS=2, nP=2) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-40,0})));

    Components.ElectronTransfer electrone
      annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                                 //(substanceData=Chemical.Examples.Substances.Electrone_solid)
    Components.ElectronTransfer electrone1
      annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                  //(substanceData=Chemical.Examples.Substances.Electrone_solid)

    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{84,-84},{104,-64}})));
    equation
      connect(Ag.port_a, electrodeReaction1.substrates[1]) annotation (Line(
          points={{-52,-20},{-42,-20},{-42,-10},{-42,-10}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Cl.port_a, electrodeReaction1.substrates[2]) annotation (Line(
          points={{-22,-26},{-38,-26},{-38,-10}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(AgCl.port_a, electrodeReaction1.products[1]) annotation (Line(
          points={{-56,14},{-42,14},{-42,10}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(H2.port_a, electrodeReaction.substrates[1]) annotation (Line(
          points={{44,42},{52,42},{52,16}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(H.port_a, electrodeReaction.products[1]) annotation (Line(
          points={{26,-26},{54,-26},{54,-4}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(electrodeReaction.products[2], electrone1.port_a) annotation (Line(
          points={{50,-4},{50,-16},{68,-16}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(electrodeReaction1.products[2], electrone.port_a) annotation (Line(
          points={{-38,10},{-38,42},{-58,42}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Cl.solution, solution1.solution) annotation (Line(
          points={{-6,-36},{-6,-40},{24.4,-40},{24.4,-59.34}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(H.solution, solution1.solution) annotation (Line(points={{10,-36},
            {10,-40},{24.4,-40},{24.4,-59.34}},
                                         smooth=Smooth.None));
    connect(electrone.solution, cathode.solution) annotation (Line(
        points={{-74,32},{-74,-42.84},{-54.4,-42.84}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(electrone1.solution, anode.solution) annotation (Line(
        points={{84,-26},{84,-49},{89.2,-49}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(AgCl.solution, cathode.solution) annotation (Line(
        points={{-72,4},{-74,4},{-74,-42.84},{-54.4,-42.84}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(Ag.solution, cathode.solution) annotation (Line(
        points={{-68,-30},{-68,-42.84},{-54.4,-42.84}},
        color={158,66,200},
        smooth=Smooth.None));
      connect(voltageSensor.p, electrone.pin) annotation (Line(
          points={{-6,74},{-96,74},{-96,42},{-78,42}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(voltageSensor.n, electrone1.pin) annotation (Line(
          points={{14,74},{92,74},{92,-16},{88,-16}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(electrone1.pin, ground.p) annotation (Line(
          points={{88,-16},{92,-16},{92,-64},{94,-64}},
          color={0,0,255},
          smooth=Smooth.None));
      annotation (
      experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ElectrochemicalCell;

    model LeadAcidBattery
    "The electrochemical cell: PbSO4(s) | Pb(s) | HSO4-(aq) , H+(aq) | PbO2(s) | PbSO4(s) + 2 H2O"
     extends Modelica.Icons.Example;

    Components.SimpleSolution anode( ElectricGround=
          false)
      annotation (Placement(transformation(extent={{24,-76},{58,32}})));

    Components.SimpleSolution cathode(ElectricGround=false)
      annotation (Placement(transformation(extent={{-80,-78},{-46,30}})));

    Components.SimpleSolution solution1(ElectricGround=false)
      annotation (Placement(transformation(extent={{-26,-80},{2,20}})));

      Components.Substance  Pb(substanceData=Chemical.Examples.Substances.Lead_solid,
        amountOfSubstance_start=1)
        annotation (Placement(transformation(extent={{50,-66},{30,-46}})));
      Components.Substance HSO4(                             substanceData=Chemical.Examples.Substances.HydrogenSulfate_aqueous,
        amountOfSubstance_start=1)
        annotation (Placement(transformation(extent={{-2,-70},{-22,-50}})));
      Components.Substance  PbSO4_(substanceData=Chemical.Examples.Substances.LeadSulfate_solid,
        amountOfSubstance_start=0.01)
        annotation (Placement(transformation(extent={{50,-32},{30,-12}})));
      Components.Substance H(substanceData=
            Chemical.Examples.Substances.Proton_aqueous,
        amountOfSubstance_start=1)    annotation (Placement(transformation(extent={{-2,-42},
              {-22,-22}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-32,72},{-12,92}})));
      Components.Reaction electrodeReaction(nP=2,
        nS=4,
        s={1,1,3,2},
      p={1,2})                                             annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-36,-14})));
      Components.Reaction electrodeReaction1(nS=2,
        nP=3,
        p={1,1,2})                                       annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={14,-16})));

    Components.ElectronTransfer electrone
      annotation (Placement(transformation(extent={{50,2},{30,22}})));
    Components.ElectronTransfer electrone1
      annotation (Placement(transformation(extent={{-72,-38},{-52,-18}})));
      Components.Substance  PbO2(substanceData=Chemical.Examples.Substances.LeadDioxide_solid,
        amountOfSubstance_start=1)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            origin={-60,-58})));
      Components.Substance H2O(substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=0.1)
        annotation (Placement(transformation(extent={{-2,-8},{-22,12}})));
      Components.Substance  PbSO4(substanceData=Chemical.Examples.Substances.LeadSulfate_solid,
        amountOfSubstance_start=0.01)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            origin={-60,6})));

    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{16,30},{36,50}})));
    Modelica.Electrical.Analog.Basic.Resistor resistor(R=1)
      annotation (Placement(transformation(extent={{-14,40},{6,60}})));
    Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
      annotation (Placement(transformation(extent={{-56,40},{-36,60}})));
    equation
      connect(Pb.port_a, electrodeReaction1.substrates[1]) annotation (Line(
          points={{30,-56},{15.5,-56},{15.5,-26},{16,-26}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(HSO4.port_a, electrodeReaction1.substrates[2]) annotation (Line(
          points={{-22,-60},{12,-60},{12,-26},{12,-26}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(PbSO4_.port_a, electrodeReaction1.products[1]) annotation (Line(
          points={{30,-22},{26,-22},{26,-2},{16,-2},{16,-6},{16.6667,-6}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(HSO4.solution, solution1.solution) annotation (Line(
          points={{-6,-70},{-6,-70},{-6,-78},{-3.6,-78},{-3.6,-79}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(H.solution, solution1.solution) annotation (Line(points={{-6,-42},
            {-6,-78},{-3.6,-78},{-3.6,-79}},
                                         smooth=Smooth.None));
      connect(H2O.solution, solution1.solution) annotation (Line(
          points={{-6,-8},{-6,-79},{-3.6,-79}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(electrodeReaction.products[1], PbSO4.port_a) annotation (Line(
          points={{-38,-4},{-38,6},{-50,6}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(electrodeReaction.products[2], H2O.port_a) annotation (Line(
          points={{-34,-4},{-34,-4},{-34,2},{-22,2}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(PbO2.port_a, electrodeReaction.substrates[1]) annotation (Line(
          points={{-50,-58},{-36,-58},{-36,-24},{-39,-24}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(HSO4.port_a, electrodeReaction.substrates[2]) annotation (Line(
          points={{-22,-60},{-34,-60},{-34,-24},{-37,-24}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(H.port_a, electrodeReaction.substrates[3]) annotation (Line(
          points={{-22,-32},{-32,-32},{-32,-24},{-35,-24}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(electrone1.port_a, electrodeReaction.substrates[4]) annotation (Line(
          points={{-52,-28},{-38,-28},{-38,-24},{-33,-24}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(
          points={{-22,-32},{2,-32},{2,2},{12,2},{12,-6},{14,-6}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
          points={{30,12},{14,12},{14,-6},{11.3333,-6}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
    connect(Pb.solution, anode.solution) annotation (Line(
        points={{46,-66},{46,-74.92},{51.2,-74.92}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(PbSO4_.solution, anode.solution) annotation (Line(
        points={{46,-32},{46,-74.92},{51.2,-74.92}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(PbO2.solution, cathode.solution) annotation (Line(
        points={{-66,-68},{-66,-76.92},{-52.8,-76.92}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(PbSO4_.solution, Pb.solution) annotation (Line(
        points={{46,-32},{46,-66}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(electrone1.pin, voltageSensor.p) annotation (Line(
        points={{-72,-28},{-82,-28},{-82,50},{-64,50},{-64,82},{-32,82}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(electrone.pin, voltageSensor.n) annotation (Line(
        points={{50,12},{50,50},{26,50},{26,82},{-12,82},{-12,82}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(electrone.solution, anode.solution) annotation (Line(
        points={{46,2},{46,-74.92},{51.2,-74.92}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(electrone.pin, ground.p) annotation (Line(
        points={{50,12},{50,50},{26,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(electrone1.pin, currentSensor.p) annotation (Line(
        points={{-72,-28},{-82,-28},{-82,50},{-56,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(currentSensor.n, resistor.p) annotation (Line(
        points={{-36,50},{-14,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(resistor.n, electrone.pin) annotation (Line(
        points={{6,50},{50,50},{50,12}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(PbSO4.solution, cathode.solution) annotation (Line(
        points={{-66,-4},{-66,-76.92},{-52.8,-76.92}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(electrone1.solution, cathode.solution) annotation (Line(
        points={{-68,-38},{-66,-38},{-66,-76.92},{-52.8,-76.92}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(H2O.solution, H.solution) annotation (Line(
        points={{-6,-8},{-6,-42}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(resistor.n, ground.p) annotation (Line(
        points={{6,50},{26,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(electrone.solution, PbSO4_.solution) annotation (Line(
        points={{46,2},{46,-32}},
        color={158,66,200},
        smooth=Smooth.None));
      annotation (
      experiment(StopTime=49500), Documentation(revisions=
                        "<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end LeadAcidBattery;

    model RedCellMembrane
     // import Chemical;
      extends Modelica.Icons.Example;

      parameter Real KC=1;//e-6 "Slow down factor";

    Chemical.Components.SimpleSolution blood_erythrocytes(
                                     ElectricGround=false)
      annotation (Placement(transformation(extent={{-180,-100},{180,-10}})));
    Chemical.Components.SimpleSolution blood_plasma
      annotation (Placement(transformation(extent={{-180,12},{180,100}})));

      Components.Substance HCO3(                               substanceData=
          Chemical.Examples.Substances.Bicarbonate_blood, amountOfSubstance_start(
            displayUnit="mmol") = 0.024)                    annotation (
        Placement(transformation(extent={{-10,-10},{10,10}}, origin={-18,30})));

      Components.Substance H2O(                            substanceData=
          Chemical.Examples.Substances.Water_liquid, amountOfSubstance_start=
          51.8*0.994648)
      annotation (Placement(transformation(extent={{-146,20},{-166,40}})));
      Components.Substance HCO3_E(      substanceData=Chemical.Examples.Substances.Bicarbonate_blood,
          amountOfSubstance_start(displayUnit="mmol") = 0.0116)
        annotation (Placement(transformation(extent={{-28,-38},{-8,-18}})));
      Components.Substance H2O_E(
        substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=38.7*0.994648)
        annotation (Placement(transformation(extent={{-144,-38},{-164,-18}})));
      Components.Substance Cl_E(
        substanceData=Chemical.Examples.Substances.Chloride_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.0499)
        annotation (Placement(transformation(extent={{-4,-38},{16,-18}})));
      Components.Substance Cl(                               substanceData=
          Chemical.Examples.Substances.Chloride_aqueous, amountOfSubstance_start(
            displayUnit="mmol") = 0.103)
      annotation (Placement(transformation(extent={{-4,20},{16,40}})));

    //  Real pH_e; //,pH_p;

      Components.Membrane Aquapirin(KC=KC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-168,0})));
      Components.Membrane Band3(KC=KC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-6,0})));
      Components.Membrane Band3_(useKineticsInput=false, KC=KC) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={18,0})));
      Components.Substance permeableUncharged(amountOfSubstance_start(displayUnit="mmol")=
             0.0118)
        annotation (Placement(transformation(extent={{166,20},{146,40}})));
      Components.Substance permeableUncharged_E(amountOfSubstance_start(displayUnit=
             "mmol") = 0.00903, substanceData(MolarWeight=0.1))
        annotation (Placement(transformation(extent={{164,-38},{144,-18}})));
      Components.Substance chargedImpermeable_E(
          amountOfSubstance_start(displayUnit="mmol") = 0.0165, substanceData(
            MolarWeight=1))
        annotation (Placement(transformation(extent={{144,-62},{164,-42}})));
      Components.Membrane leak(useKineticsInput=false, KC=KC) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={140,0})));
      Components.Substance Lac_E(
        substanceData=Chemical.Examples.Substances.Chloride_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.00062)
        annotation (Placement(transformation(extent={{56,-38},{76,-18}})));
      Components.Substance Lac(substanceData=Chemical.Examples.Substances.Chloride_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.00131)
        annotation (Placement(transformation(extent={{56,20},{76,40}})));
      Components.Membrane MCT_(useKineticsInput=false, KC=KC)
      "Monocarboxylate transporters"   annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={78,0})));
      Components.Substance H_E(substanceData=Chemical.Examples.Substances.Proton_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 38.7*10^(-7.2)) "H+"
        annotation (Placement(transformation(extent={{30,-38},{50,-18}})));
      Components.Substance H(substanceData=Chemical.Examples.Substances.Proton_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 51.8*10^(-7.4))
      "H+ in plasma"
        annotation (Placement(transformation(extent={{30,20},{50,40}})));
      Components.Membrane MCT(useKineticsInput=false, KC=KC)
      "Monocarboxylate transporters"   annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={52,0})));
      Components.Substance CO2(substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.00167)
      "free dissolved unbound CO2"
      annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
      Components.Substance CO2_E(substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.00125)
      "free dissolved unbound CO2"
        annotation (Placement(transformation(extent={{-58,-38},{-38,-18}})));
      Components.Membrane freeCO2(KC=KC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-36,0})));
      Components.Substance O2(substanceData=Chemical.Examples.Substances.Oxygen_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.000167)
      "free dissolved undound oxygen"
        annotation (Placement(transformation(extent={{96,20},{116,40}})));
      Components.Membrane freeO2(KC=KC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={118,0})));
      Components.Substance O2_E(amountOfSubstance_start(displayUnit="mmol") = 0.000125,
          substanceData=Chemical.Examples.Substances.Oxygen_aqueous)
      "free dissolved undound O2"
        annotation (Placement(transformation(extent={{96,-38},{116,-18}})));
      Chemical.Components.Substance
                           K(substanceData=Chemical.Examples.Substances.Potassium_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.004)
        annotation (Placement(transformation(extent={{-92,20},{-112,40}})));
      Chemical.Components.Substance
                           Na(                               substanceData=Chemical.Examples.Substances.Sodium_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.138)
        annotation (Placement(transformation(extent={{-118,20},{-138,40}})));
      Chemical.Components.Substance
                           Na_E(substanceData=Chemical.Examples.Substances.Sodium_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.007)
        annotation (Placement(transformation(extent={{-118,-38},{-138,-18}})));
      Chemical.Components.Substance
                           K_E(substanceData=Chemical.Examples.Substances.Potassium_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.096)
        annotation (Placement(transformation(extent={{-112,-38},{-92,-18}})));
      Chemical.Components.Substance H2PO4_E(substanceData=Chemical.Examples.Substances.DihydrogenPhosphate_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.000175)
        annotation (Placement(transformation(extent={{-84,-38},{-64,-18}})));
      Chemical.Components.Substance ADP_E(substanceData(z=-3),
          amountOfSubstance_start(displayUnit="mmol") = 9.6e-05)
        annotation (Placement(transformation(extent={{-114,-62},{-94,-42}})));
      Chemical.Components.Substance ATP_E(substanceData(
          z=-4,
          DfH_25degC=16700,
          DfG_25degC_1bar=30500,
          References={"http://www.wiley.com/college/pratt/0471393878/student/review/thermodynamics/7_relationship.html"}),
          amountOfSubstance_start(displayUnit="mmol") = 0.00128)
        annotation (Placement(transformation(extent={{-118,-62},{-138,-42}})));
      Chemical.Components.Substance HPO4_E(substanceData=Chemical.Examples.Substances.HydrogenPhosphate_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.000495)
        annotation (Placement(transformation(extent={{-84,-62},{-64,-42}})));
      Chemical.Components.Substance H2PO4(substanceData=Chemical.Examples.Substances.DihydrogenPhosphate_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.000365)
        annotation (Placement(transformation(extent={{-86,20},{-66,40}})));
      Chemical.Components.Substance HPO4(substanceData=Chemical.Examples.Substances.HydrogenPhosphate_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.001635)
        annotation (Placement(transformation(extent={{-86,42},{-66,62}})));
      Components.Substance albumin(substanceData(
          MolarWeight=66.463,
          z=-17,
          density=1080), amountOfSubstance_start(displayUnit="mmol") = 0.0007)
        annotation (Placement(transformation(extent={{116,76},{96,96}})));
      Components.Substance globulins(substanceData(
          MolarWeight=34,
          z=-2.43,
          density=1080), amountOfSubstance_start(displayUnit="mmol") = 0.00082)
        annotation (Placement(transformation(extent={{150,76},{130,96}})));
      Chemical.Components.Substance Ca(substanceData=Chemical.Examples.Substances.Calcium_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.00175) "Ca2+"
        annotation (Placement(transformation(extent={{-112,42},{-92,62}})));
      Chemical.Components.Substance Mg(substanceData=Chemical.Examples.Substances.Magnesium_aqueous,
          amountOfSubstance_start(displayUnit="mmol") = 0.00108) "Mg2+"
        annotation (Placement(transformation(extent={{-112,-84},{-92,-64}})));
      Components.Substance hemoglobin(substanceData(
          MolarWeight=64,
          z=-4.4,
          density=1500), amountOfSubstance_start(displayUnit="mmol") = 0.00513)
        annotation (Placement(transformation(extent={{94,-94},{74,-74}})));
      Components.Substance DPG(amountOfSubstance_start(displayUnit="mmol") = 0.0051,
          substanceData(
          MolarWeight=0.266,
          z=-2.2,
          density=1000))
        annotation (Placement(transformation(extent={{128,-94},{108,-74}})));
      Components.Substance GSH(substanceData(
          MolarWeight=0.2,
          z=-1,
          density=1000), amountOfSubstance_start(displayUnit="mmol") = 0.00223)
        annotation (Placement(transformation(extent={{164,-94},{144,-74}})));
    equation
    //  pH_p = -log10(H.a);
     // pH_e = -log10(H_E.a);
    connect(H2O.solution, blood_plasma.solution)
      annotation (Line(points={{-150,20},{108,20},{108,12.88}},
                                                        smooth=Smooth.None));
    connect(Cl.solution, blood_plasma.solution) annotation (Line(
        points={{0,20},{0,16},{0,12.88},{108,12.88}},
        color={0,0,0},
        smooth=Smooth.None));
      connect(H2O_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{-148,-38},{108,-38},{108,-99.1}},
                                                    smooth=Smooth.None));
      connect(Cl_E.solution, blood_erythrocytes.solution) annotation (Line(
          points={{0,-38},{0,-68},{0,-99.1},{108,-99.1}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(HCO3_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{-24,-38},{108,-38},{108,-99.1}},
                                                  smooth=Smooth.None));
      connect(Aquapirin.port_b, H2O_E.port_a) annotation (Line(
          points={{-168,-10},{-168,-28},{-164,-28}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Aquapirin.port_a, H2O.port_a) annotation (Line(
          points={{-168,10},{-168,30},{-166,30}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Band3.port_a, HCO3.port_a) annotation (Line(
          points={{-6,10},{-6,30},{-8,30}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Band3.port_b, HCO3_E.port_a) annotation (Line(
          points={{-6,-10},{-6,-28},{-8,-28}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Band3_.port_b, Cl_E.port_a) annotation (Line(
          points={{18,-10},{18,-28},{16,-28}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Band3_.port_a, Cl.port_a) annotation (Line(
          points={{18,10},{18,30},{16,30}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
    connect(HCO3.solution, blood_plasma.solution) annotation (Line(
        points={{-24,20},{108,20},{108,12.88}},
        color={0,0,0},
        smooth=Smooth.None));
      connect(blood_plasma.solution, permeableUncharged.solution) annotation (Line(
          points={{108,12.88},{108,20},{162,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(blood_erythrocytes.solution, permeableUncharged_E.solution)
        annotation (Line(
          points={{108,-99.1},{108,-38},{160,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(blood_erythrocytes.solution,chargedImpermeable_E. solution)
        annotation (Line(
          points={{108,-99.1},{108,-38},{140,-38},{140,-62},{148,-62}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(permeableUncharged.port_a, leak.port_a) annotation (Line(
          points={{146,30},{140,30},{140,10}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(permeableUncharged_E.port_a, leak.port_b) annotation (Line(
          points={{144,-28},{140,-28},{140,-10}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(MCT_.port_a, Lac.port_a) annotation (Line(
          points={{78,10},{78,30},{76,30}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(MCT_.port_b, Lac_E.port_a) annotation (Line(
          points={{78,-10},{78,-28},{76,-28}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(Lac.solution, blood_plasma.solution) annotation (Line(
          points={{60,20},{108,20},{108,12.88}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(blood_erythrocytes.solution, Lac_E.solution) annotation (Line(
          points={{108,-99.1},{108,-38},{60,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H_E.solution, blood_erythrocytes.solution) annotation (Line(
          points={{34,-38},{108,-38},{108,-99.1}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H_E.port_a, MCT.port_b) annotation (Line(
          points={{50,-28},{52,-28},{52,-10}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(MCT.port_a, H.port_a) annotation (Line(
          points={{52,10},{52,30},{50,30}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(blood_plasma.solution, H.solution) annotation (Line(
          points={{108,12.88},{108,20},{34,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(CO2.port_a, freeCO2.port_a) annotation (Line(
          points={{-40,30},{-36,30},{-36,10}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(freeCO2.port_b, CO2_E.port_a) annotation (Line(
          points={{-36,-10},{-36,-28},{-38,-28}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(blood_plasma.solution, CO2.solution) annotation (Line(
          points={{108,12.88},{108,20},{-56,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(CO2_E.solution, blood_erythrocytes.solution) annotation (Line(
          points={{-54,-38},{108,-38},{108,-99.1}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(blood_plasma.solution, O2.solution) annotation (Line(
          points={{108,12.88},{108,20},{100,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(O2_E.solution, blood_erythrocytes.solution) annotation (Line(
          points={{100,-38},{108,-38},{108,-99.1}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(O2_E.port_a, freeO2.port_b) annotation (Line(
          points={{116,-28},{118,-28},{118,-10}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(freeO2.port_a, O2.port_a) annotation (Line(
          points={{118,10},{118,30},{116,30}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(H2O.solution, K.solution) annotation (Line(
          points={{-150,20},{-96,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H2O.solution, Na.solution) annotation (Line(
          points={{-150,20},{-122,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H2O_E.solution, Na_E.solution) annotation (Line(
          points={{-148,-38},{-122,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H2O_E.solution, K_E.solution) annotation (Line(
          points={{-148,-38},{-108,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H2O_E.solution, H2PO4_E.solution) annotation (Line(
          points={{-148,-38},{-80,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(ADP_E.solution, K_E.solution) annotation (Line(
          points={{-110,-62},{-110,-38},{-108,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(ATP_E.solution, Na_E.solution) annotation (Line(
          points={{-122,-62},{-122,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H2O.solution, H2PO4.solution) annotation (Line(
          points={{-150,20},{-82,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(HPO4_E.solution, H2PO4_E.solution) annotation (Line(
          points={{-80,-62},{-110,-62},{-110,-38},{-80,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(HPO4.solution, H2PO4.solution) annotation (Line(
          points={{-82,42},{-82,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(albumin.solution, permeableUncharged.solution) annotation (Line(
          points={{112,76},{92,76},{92,20},{162,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(globulins.solution, permeableUncharged.solution) annotation (Line(
          points={{146,76},{92,76},{92,20},{162,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(Ca.solution, CO2.solution) annotation (Line(
          points={{-108,42},{-82,42},{-82,20},{-56,20}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(Mg.solution, blood_erythrocytes.solution) annotation (Line(
          points={{-108,-84},{-108,-38},{108,-38},{108,-99.1}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(DPG.solution, permeableUncharged_E.solution) annotation (Line(
          points={{124,-94},{140,-94},{140,-38},{160,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(hemoglobin.solution, permeableUncharged_E.solution) annotation (Line(
          points={{90,-94},{140,-94},{140,-38},{160,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(GSH.solution, permeableUncharged_E.solution) annotation (Line(
          points={{160,-94},{140,-94},{140,-38},{160,-38}},
          color={158,66,200},
          smooth=Smooth.None));
      annotation ( Documentation(info="<html>
<p>Blood eqiulibrium across erythrocyte membrane bewteen blood plasma and intracellular fluid of erythrocytes.</p>
<p>Data of blood status are from:</p>
<p>Raftos, J.E., Bulliman, B.T. and Kuchel, P.W. Evaluation of an electrochemical model of erythrocyte pH buffering using 31P nuclear magnetic resonance data. <i>The Journal of general physiology</i> 1990;95(6):1183-1204. </p>
</html>",    revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment(StopTime=1e-008),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},{180,100}}),
                           graphics));
    end RedCellMembrane;

    package AcidBase

      model WaterSelfIonization "H2O  <->  OH-   +   H+ "
        import Chemical;
          extends Modelica.Icons.Example;

      Chemical.Components.SimpleSolution solution
        annotation (Placement(transformation(extent={{-72,2},{76,96}})));
      Chemical.Components.SimpleSolution solution1
        annotation (Placement(transformation(extent={{-74,-96},{74,-2}})));
        Components.Substance H3O(
          amountOfSubstance_start=1e-7,   substanceData=
              Chemical.Examples.Substances.Hydronium_aqueous)
                                        annotation (Placement(transformation(extent={{10,-10},
                {-10,10}},       origin={30,70})));
        Components.Substance OH(
          amountOfSubstance_start=1e-7,   substanceData=
              Chemical.Examples.Substances.Hydroxide_aqueous)
                                        annotation (Placement(transformation(extent={{10,-10},
                {-10,10}},       origin={30,26})));
        Components.Substance H2O(
          amountOfSubstance_start=1,   substanceData=
              Chemical.Examples.Substances.Water_liquid)
                                     annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={-30,46})));
        Chemical.Components.Reaction waterDissociation(nP=2, s={2})
          annotation (Placement(transformation(extent={{-12,36},{8,56}})));
              Real pH, pH_;
        Components.Substance H_(
          amountOfSubstance_start=1e-7,   substanceData=
              Chemical.Examples.Substances.Proton_aqueous)       annotation (Placement(transformation(extent={{10,-10},
                {-10,10}},       origin={28,-30})));
        Components.Substance OH_(
          amountOfSubstance_start=1e-7,   substanceData=
              Chemical.Examples.Substances.Hydroxide_aqueous)
                                        annotation (Placement(transformation(extent={{10,-10},
                {-10,10}},       origin={28,-76})));
        Components.Substance H2O_(
          amountOfSubstance_start=1,   substanceData=
              Chemical.Examples.Substances.Water_liquid)
                                     annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={-32,-56})));
        Chemical.Components.Reaction waterDissociation_(nP=2)
          annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));

      equation
        pH = -log10( H3O.a);

        pH_ = -log10( H_.a);

        connect(OH.port_a, waterDissociation.products[1]) annotation (Line(
            points={{20,26},{16,26},{16,44},{8,44}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(waterDissociation.products[2], H3O.port_a) annotation (Line(
            points={{8,48},{16,48},{16,70},{20,70}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.port_a, waterDissociation.substrates[1]) annotation (Line(
            points={{-20,46},{-12,46}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(OH_.port_a,waterDissociation_. products[1]) annotation (Line(
            points={{18,-76},{14,-76},{14,-58},{6,-58}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(waterDissociation_.products[2], H_.port_a) annotation (Line(
            points={{6,-54},{14,-54},{14,-30},{18,-30}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O_.port_a,waterDissociation_. substrates[1]) annotation (Line(
            points={{-22,-56},{-14,-56}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.solution, solution.solution) annotation (Line(
            points={{-36,36},{2,36},{2,2.94}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(OH.solution, solution.solution) annotation (Line(
            points={{36,16},{36,2.94},{2,2.94}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H3O.solution, solution.solution) annotation (Line(
            points={{36,60},{36,2.94},{2,2.94}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H2O_.solution, solution1.solution) annotation (Line(
            points={{-38,-66},{0,-66},{0,-95.06}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(OH_.solution, solution1.solution) annotation (Line(
            points={{34,-86},{34,-95.06},{0,-95.06}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H_.solution, solution1.solution) annotation (Line(
            points={{34,-40},{34,-95.06},{0,-95.06}},
            color={0,0,0},
            smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>Self-ionization of water.</p>
<p>Ions difference (SID) in water causes the acidity/basicity, where pH = -log10(aH+). An activity of hydrogen ions aH+ is approximated with concentration (mol/l) of the oxonium cations H3O+.</p>
<pre><b>plotExpression(apply(-log10(WaterSelfIonization.H3O.solute)),&nbsp;false,&nbsp;&QUOT;pH&QUOT;,&nbsp;1);</b></pre>
<p><br>The titration slope der(pH)/der(SID)=1.48e+6 1/(mol/L) at pH=7.4.</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=1));
      end WaterSelfIonization;

      model CarbonDioxideInWater "CO2 as alone acid-base buffer"
        import Chemical;
          extends Modelica.Icons.Example;
      Chemical.Components.SimpleSolution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,46}})));
        Components.Substance HCO3(  substanceData=
              Chemical.Examples.Substances.Bicarbonate_aqueous)
          annotation (Placement(transformation(extent={{-16,-4},{4,16}})));
        Chemical.Components.Reaction HendersonHasselbalch(nP=2, nS=2,
        useKineticsInput=false) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-50,-6},{-30,14}})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
        substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
        PartialPressure=5332.8954966,
        TotalPressure=101325.0144354) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-60,86})));
        Components.Substance H(  substanceData=
              Chemical.Examples.Substances.Proton_aqueous,
            amountOfSubstance_start=3e-8) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}}, origin={-6,-38})));
        Components.GasSolubility gasSolubility
          annotation (Placement(transformation(extent={{-70,36},{-50,56}})));
                                              /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
        Components.Substance CO2_liquid(  substanceData=
              Chemical.Examples.Substances.CarbonDioxide_aqueous)
          annotation (Placement(transformation(extent={{-82,-6},{-62,14}})));
        Components.Substance CO3(  substanceData=
              Chemical.Examples.Substances.Carbonate_aqueous)
          annotation (Placement(transformation(extent={{70,-2},{50,18}})));
        Chemical.Components.Reaction c2(nP=2, nS=1)
        "K=10^(-10.33 + 3), dH=14.9kJ/mol"
          annotation (Placement(transformation(extent={{16,-4},{36,16}})));
        Chemical.Components.Substance H2O(  substanceData=
              Chemical.Examples.Substances.Water_liquid,
            amountOfSubstance_start=55.507)
          annotation (Placement(transformation(extent={{-82,-50},{-62,-30}})));
        Real pH;

      equation
        pH = -log10( H.a);

        connect(CO2_gas.port_a, gasSolubility.gas_port) annotation (Line(
            points={{-60,76},{-60,56}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(gasSolubility.liquid_port, CO2_liquid.port_a) annotation (Line(
            points={{-60,36},{-60,4},{-62,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch.products[1], H.port_a) annotation (Line(
            points={{-30,2},{-22,2},{-22,-38},{4,-38}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HendersonHasselbalch.products[2], HCO3.port_a) annotation (Line(
            points={{-30,6},{-22,6},{-22,6},{4,6}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HCO3.port_a, c2.substrates[1]) annotation (Line(
            points={{4,6},{16,6}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(c2.products[1], H.port_a) annotation (Line(
            points={{36,4},{44,4},{44,-38},{4,-38}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(c2.products[2], CO3.port_a) annotation (Line(
            points={{36,8},{48,8},{48,8},{50,8}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.port_a, HendersonHasselbalch.substrates[2]) annotation (
            Line(
            points={{-62,4},{-62,6},{-50,6}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2O.port_a, HendersonHasselbalch.substrates[1]) annotation (Line(
            points={{-62,-40},{-56,-40},{-56,2},{-50,2}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(CO2_liquid.solution, solution.solution) annotation (Line(
            points={{-78,-6},{-78,-98.54},{0,-98.54}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H2O.solution, solution.solution) annotation (Line(points={{-78,-50},
              {-78,-98.54},{0,-98.54}}, smooth=Smooth.None));
        connect(HCO3.solution, solution.solution) annotation (Line(points={{-12,-4},
              {-12,-98.54},{0,-98.54}}, smooth=Smooth.None));
        connect(H.solution, solution.solution) annotation (Line(points={{-12,-48},
              {-12,-98.54},{0,-98.54}}, smooth=Smooth.None));
        connect(CO3.solution, solution.solution) annotation (Line(points={{66,-2},
              {66,-98.54},{0,-98.54}},  smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>CO2 solution in water without any other acid-base buffers.</p>
<pre><b>plotExpression(apply(-log10(CarbonDioxideInWater.H3O.solute)),&nbsp;false,&nbsp;&QUOT;pH&QUOT;,&nbsp;1);</b></pre>
<p><br>Please note, that OH- (and CO3^-2) can be neglected from electroneutrality calculation, because of very small concentrations (in physiological pH) anyway. </p>
<p>And if SID&GT;0 then also H3O+ can be also neglected from electroneutrality, because only bicarbonate anions HCO3- (or CO3^-2) are needed there to balance the electroneutrality.</p>
<p><br>The partial pressure of CO2 in gas are input parameter. Outputs are an amount of free dissolved CO2 in liquid and an amount of HCO3-.</p>
<p><br>The titration slope der(pH)/der(SID)=17.5 1/(mol/L) at pH=7.4 and pCO2=40 mmHg.</p>
<p><br>Molar heat of formation (aqueous):</p>
<p>CO2:        -413.5 kJ/mol  (gas: -393.5 kJ/mol )</p>
<p>H2O:        -285.8 kJ/mol</p>
<p>HCO3-:        -692.0 kJ/mol</p>
<p>CO3^-2:        -677.1 kJ/mol</p>
<p><br>Enthalphy of reaction H2O + CO2 &LT;-&GT; HCO3- + H+  :         7.3 kJ/mol</p>
<p>Enthalphy of reaction HCO3- &LT;-&GT; CO3^-2 + H+  :        14.9 kJ/mol</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.02));
      end CarbonDioxideInWater;

      model Phosphate
        import Chemical;
          extends Modelica.Icons.Example;

      Chemical.Components.SimpleSolution solution
        annotation (Placement(transformation(extent={{-98,-100},{100,100}})));

         Components.Substance H(
          amountOfSubstance_start=55.6*10^(-7.4),
            substanceData=Chemical.Examples.Substances.Proton_aqueous)
        "hydrogen ions activity"   annotation (Placement(transformation(extent=
                  {{-10,-10},{10,10}}, origin={28,-14})));

        Components.Substance H3PO4(
          amountOfSubstance_start=1e-08,
            substanceData=
              Chemical.Examples.Substances.PhosphoricAcid_aqueous)
          annotation (Placement(transformation(extent={{-90,-58},{-70,-38}})));
        Components.Substance H2PO4(
          amountOfSubstance_start=0.0005,
            substanceData=
              Chemical.Examples.Substances.DihydrogenPhosphate_aqueous)
          annotation (Placement(transformation(extent={{-40,-58},{-20,-38}})));
        Components.Substance HPO4(
            substanceData=
              Chemical.Examples.Substances.HydrogenPhosphate_aqueous,
          amountOfSubstance_start=0.0006)
          annotation (Placement(transformation(extent={{16,-58},{36,-38}})));
        Components.Substance PO4(
            substanceData=Chemical.Examples.Substances.Phosphate_aqueous,
          amountOfSubstance_start=1e-08)
          annotation (Placement(transformation(extent={{92,-58},{72,-38}})));

        Chemical.Components.Reaction chemicalReaction(nP=2) "10^(-1.915 + 3)"
          annotation (Placement(transformation(extent={{-66,-58},{-46,-38}})));
        Chemical.Components.Reaction chemicalReaction1(nP=2) "10^(-6.66 + 3)"
          annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
        Chemical.Components.Reaction chemicalReaction2(nP=2) "10^(-11.78 + 3)"
          annotation (Placement(transformation(extent={{44,-58},{64,-38}})));

        Chemical.Components.Substance
                             H2O(      substanceData=
              Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=55.508)
                                     annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={58,-76})));
      equation
        connect(H3PO4.port_a, chemicalReaction.substrates[1]) annotation (Line(
            points={{-70,-48},{-66,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction.products[1], H2PO4.port_a) annotation (Line(
            points={{-46,-50},{-42,-50},{-42,-48},{-20,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H2PO4.port_a, chemicalReaction1.substrates[1]) annotation (Line(
            points={{-20,-48},{-14,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[1], HPO4.port_a) annotation (Line(
            points={{6,-50},{16,-50},{16,-48},{36,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(HPO4.port_a, chemicalReaction2.substrates[1]) annotation (Line(
            points={{36,-48},{44,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction2.products[1], PO4.port_a) annotation (Line(
            points={{64,-50},{74,-50},{74,-48},{72,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction.products[2], H.port_a) annotation (Line(
            points={{-46,-46},{-44,-46},{-44,-32},{38,-32},{38,-14}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction1.products[2], H.port_a) annotation (Line(
            points={{6,-46},{14,-46},{14,-32},{38,-32},{38,-14}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(chemicalReaction2.products[2], H.port_a) annotation (Line(
            points={{64,-46},{66,-46},{66,-32},{38,-32},{38,-14}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(H3PO4.solution, solution.solution) annotation (Line(
            points={{-86,-58},{-46,-58},{-46,-98},{1,-98}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(H2PO4.solution, solution.solution) annotation (Line(points={{-36,-58},
              {-36,-88},{1,-88},{1,-98}},    smooth=Smooth.None));
        connect(HPO4.solution, solution.solution) annotation (Line(points={{20,-58},
              {22,-58},{22,-88},{1,-88},{1,-98}},smooth=Smooth.None));
        connect(PO4.solution, solution.solution) annotation (Line(points={{88,-58},
              {88,-88},{1,-88},{1,-98}},smooth=Smooth.None));
        connect(H.solution, solution.solution) annotation (Line(points={{22,-24},
              {22,-88},{1,-88},{1,-98}},
                                   smooth=Smooth.None));
      connect(chemicalReaction.substrates[1], H3PO4.port_a) annotation (Line(
          points={{-66,-48},{-70,-48}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(H2O.solution, solution.solution) annotation (Line(
          points={{52,-86},{52,-98},{1,-98}},
          color={158,66,200},
          smooth=Smooth.None));
        annotation ( Documentation(info="<html>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.05));
      end Phosphate;

      model AlbuminTitration "Figge-Fencl model (22. Dec. 2007)"
        extends Modelica.Icons.Example;

      Components.SimpleSolution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Sources.Buffer H(substanceData=Chemical.Examples.Substances.Proton_aqueous,
          a_start=10^(-7.4)) "hydrogen ions activity" annotation (Placement(
            transformation(extent={{10,-10},{-10,10}}, origin={14,36})));

        constant Integer n=218 "Number of weak acid group in albumin molecule";
        constant Real pKAs[n]=cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
        "acid dissociation constants";
        constant Real K[n]=fill(10.0, n) .^ (-pKAs);
        constant Real DfG[n]= Modelica.Constants.R*(298.15)*log(K);

        Chemical.Components.Substance A[n](
          each amountOfSubstance_start=0.00033, substanceData(each z=-1))
        "deprotonated acid groups"
          annotation (Placement(transformation(extent={{24,-16},{4,4}})));
        Chemical.Components.Reaction react[n](each nP=2)
          annotation (Placement(transformation(extent={{-44,-2},{-24,18}})));

        Chemical.Components.Substance HA[n](substanceData(DfG_25degC_1bar=DfG), each amountOfSubstance_start=
             0.00033) "protonated acid groups"
          annotation (Placement(transformation(extent={{-78,-2},{-58,18}})));

        Components.Substance H2O(      substanceData=
              Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=55.508)
                                     annotation (Placement(transformation(extent={{-10,
                  -10},{10,10}}, origin={62,-68})));
      equation
        connect(react.products[1], A.port_a) annotation (Line(
            points={{-24,6},{-12,6},{-12,-6},{4,-6}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        for i in 1:n loop
          connect(react[i].products[2], H.port_a) annotation (Line(
              points={{-24,10},{-14,10},{-14,36},{4,36}},
              color={107,45,134},
              thickness=1,
              smooth=Smooth.None));
          connect(HA[i].solution, solution.solution) annotation (Line(
            points={{-74,-2},{-74,-86},{0,-86},{0,-98}},
            color={0,0,0},
            smooth=Smooth.None));
          connect(A[i].solution, solution.solution) annotation (Line(
            points={{20,-16},{20,-86},{0,-86},{0,-98}},
            color={0,0,0},
            smooth=Smooth.None));
        end for;
        connect(HA.port_a, react.substrates[1]) annotation (Line(
            points={{-58,8},{-44,8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

      connect(solution.solution, H2O.solution) annotation (Line(
          points={{0,-98},{56,-98},{56,-78}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H.solution, solution.solution) annotation (Line(
          points={{20,26},{20,14},{36,14},{36,-98},{0,-98}},
          color={158,66,200},
          smooth=Smooth.None));
        annotation ( Documentation(revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",       info="<html>
<p>The titration slope der(pH)/der(SID)=185 1/(mol/L) at pH=7.4 and tAlb=0.66 mmol/l.</p>
<p>Data and model is described in</p>
<p><font style=\"color: #222222; \">Jame Figge: Role of non-volatile weak acids (albumin, phosphate and citrate). In: Stewart&apos;s Textbook of Acid-Base, 2nd Edition, John A. Kellum, Paul WG Elbers editors, &nbsp;AcidBase org, 2009, pp. 216-232.</font></p>
</html>"),experiment(StopTime=1.6));
      end AlbuminTitration;

      model CarbonDioxideInBlood
        import Chemical;
          extends Modelica.Icons.Example;

        parameter Real KC=10;//e-6 "Slow down factor";

      Chemical.Components.SimpleSolution blood_erythrocytes(ElectricGround=
            false, temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,-96},{100,-36}})));
      Chemical.Components.SimpleSolution blood_plasma(temperature_start=310.15)
        annotation (Placement(transformation(extent={{-100,6},{100,58}})));

        Components.Substance HCO3(amountOfSubstance_start=0.024, substanceData=
            Chemical.Examples.Substances.Bicarbonate_blood)   annotation (
          Placement(transformation(extent={{10,-10},{-10,10}}, origin={18,24})));
      Chemical.Sources.ExternalIdealGasSubstance CO2_gas(
        substanceData=Chemical.Examples.Substances.CarbonDioxide_gas,
        TotalPressure(displayUnit="mmHg") = 101325.0144354,
        PartialPressure(displayUnit="mmHg") = 5332.8954966,
        usePartialPressureInput=true,
        Temperature=310.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-84,84})));
        Components.GasSolubility gasSolubility(KC=KC)
          annotation (Placement(transformation(extent={{-94,48},{-74,68}})));

        Components.Substance CO2(substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous,
          amountOfSubstance_start=0.0017) "Free dissolved CO2 in plasma"
        annotation (Placement(transformation(extent={{-88,28},{-68,48}})));
        Components.Substance H2O(                            substanceData=
            Chemical.Examples.Substances.Water_liquid, amountOfSubstance_start=
            51.8*0.994648)
        annotation (Placement(transformation(extent={{-60,12},{-40,32}})));
        Components.Substance HCO3_E(
          amountOfSubstance_start=0.0116, substanceData=Chemical.Examples.Substances.Bicarbonate_blood)
          annotation (Placement(transformation(extent={{28,-62},{8,-42}})));
        Chemical.Components.Reaction HendersonHasselbalch1(nP=2, nS=2,
        KC=KC) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-26,-68},{-6,-48}})));
        Components.Substance CO2_E(substanceData=Chemical.Examples.Substances.CarbonDioxide_aqueous,
          amountOfSubstance_start=0.00123) "Free dissolved CO2 in erythrocyte"
        annotation (Placement(transformation(extent={{-90,-82},{-70,-62}})));
        Components.Substance H2O_E(
          substanceData=Chemical.Examples.Substances.Water_liquid,
            amountOfSubstance_start=38.7*0.994648)
          annotation (Placement(transformation(extent={{-60,-62},{-40,-42}})));
        Components.Substance Cl_E(
          amountOfSubstance_start=0.0499,
          substanceData=Chemical.Examples.Substances.Chloride_aqueous)
          annotation (Placement(transformation(extent={{66,-60},{46,-40}})));
        Components.Substance Cl(amountOfSubstance_start=0.103, substanceData=
            Chemical.Examples.Substances.Chloride_aqueous)
        annotation (Placement(transformation(extent={{66,20},{46,40}})));

        Real pH_e, pH_p;

        Components.Membrane aquaporin(KC=KC)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-34,-16})));
        Components.Membrane Band3_HCO3(KC=KC) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={4,-16})));
        Components.Membrane Band3_Cl(useKineticsInput=false, KC=KC) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={44,-16})));
      Chemical.Sources.Buffer H_E(
        substanceData=Chemical.Examples.Substances.Proton_aqueous,
        a_start=10^(-7.2),
        BufferValue=0.063)
        annotation (Placement(transformation(extent={{48,-84},{30,-66}})));
        Modelica.Blocks.Sources.Clock clock(offset=5000)
          annotation (Placement(transformation(extent={{-54,62},{-34,82}})));
        Components.Substance others_E(amountOfSubstance_start=38.7*(1 -
            0.994648) - 0.0499 - 0.0116 - 0.00123, substanceData(
          density=(1.045 - 0.695523)*1000/(1 - 0.697583),
          References={"erythrocyte intracellular fluid density 1045kg/m3"},
          MolarWeight=(1.045 - 0.695523)/(38.7*(1 - 0.994648) - 0.0499 - 0.0116
               - 0.00123)))
        annotation (Placement(transformation(extent={{68,-88},{88,-68}})));
        Components.Substance others_P(amountOfSubstance_start=51.8*(1 -
            0.994648) - 0.103 - 0.024 - 0.0017, substanceData(
          References={
              "to reach plasma density 1024 kg/m3 and plasma volume 1 liter"},
          density=(1.024 - 0.933373)*1000/(1 - 0.936137),
          MolarWeight=(1.024 - 0.933373)/(51.8*(1 - 0.994648) - 0.103 - 0.024
               - 0.0017)))
        annotation (Placement(transformation(extent={{70,14},{90,34}})));
        Chemical.Components.Diffusion diffusion annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-66,-16})));
      Chemical.Sources.Buffer H(
        substanceData=Chemical.Examples.Substances.Proton_aqueous,
        a_start=10^(-7.4),
        BufferValue=0.0077)
        "buffer value 7.7 mmol/L for plasma is from (O. Siggaard-Andersen 1995)"
        annotation (Placement(transformation(extent={{40,38},{22,56}})));
        Chemical.Components.Reaction HendersonHasselbalch2(nP=2, nS=2,
          KC=(1e-10)*KC) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-26,26},{-6,46}})));
      equation
        pH_p = -log10(H.a);
        pH_e = -log10(H_E.a);
        connect(HendersonHasselbalch1.products[1], HCO3_E.port_a) annotation (Line(
            points={{-6,-60},{2,-60},{2,-52},{8,-52}},
            color={107,45,134},
            thickness=0.5,
            smooth=Smooth.None));
      connect(CO2_E.port_a, HendersonHasselbalch1.substrates[1]) annotation (
          Line(
          points={{-70,-72},{-36,-72},{-36,-60},{-26,-60}},
          color={107,45,134},
          thickness=0.5,
          smooth=Smooth.None));
        connect(H2O_E.port_a, HendersonHasselbalch1.substrates[2]) annotation (Line(
            points={{-40,-52},{-34,-52},{-34,-56},{-26,-56}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
      connect(CO2.solution, blood_plasma.solution) annotation (Line(
          points={{-84,28},{-84,6.52},{0,6.52}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(H2O.solution, blood_plasma.solution)
        annotation (Line(points={{-56,12},{-56,6.52},{0,6.52}},
                                                          smooth=Smooth.None));
      connect(Cl.solution, blood_plasma.solution) annotation (Line(
          points={{62,20},{62,6.52},{0,6.52}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(CO2_E.solution, blood_erythrocytes.solution) annotation (Line(
          points={{-86,-82},{-86,-88},{0,-88},{0,-95.4}},
          color={0,0,0},
          smooth=Smooth.None));
        connect(H2O_E.solution, blood_erythrocytes.solution) annotation (Line(
              points={{-56,-62},{-56,-88},{0,-88},{0,-95.4}},
                                                      smooth=Smooth.None));
        connect(Cl_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{62,-60},{62,-88},{0,-88},{0,-95.4}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(HCO3_E.solution, blood_erythrocytes.solution) annotation (Line(
              points={{24,-62},{24,-88},{0,-88},{0,-95.4}},
                                                    smooth=Smooth.None));
      connect(gasSolubility.liquid_port, CO2.port_a) annotation (Line(
          points={{-84,48},{-84,38},{-68,38}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
        connect(aquaporin.port_b, H2O_E.port_a) annotation (Line(
            points={{-34,-26},{-34,-52},{-40,-52}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
      connect(aquaporin.port_a, H2O.port_a) annotation (Line(
          points={{-34,-6},{-34,22},{-40,22}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
        connect(Band3_HCO3.port_a, HCO3.port_a) annotation (Line(
            points={{4,-6},{4,24},{8,24}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        connect(Band3_HCO3.port_b, HCO3_E.port_a) annotation (Line(
            points={{4,-26},{4,-52},{8,-52}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        connect(Band3_Cl.port_b, Cl_E.port_a) annotation (Line(
            points={{44,-26},{44,-50},{46,-50}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        connect(Band3_Cl.port_a, Cl.port_a) annotation (Line(
            points={{44,-6},{44,30},{46,30}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        connect(gasSolubility.gas_port, CO2_gas.port_a) annotation (Line(
            points={{-84,68},{-84,74}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
      connect(HCO3.solution, blood_plasma.solution) annotation (Line(
          points={{24,14},{24,6.52},{0,6.52}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(H_E.port_a, HendersonHasselbalch1.products[2]) annotation (Line(
          points={{30,-75},{4,-75},{4,-56},{-6,-56}},
          color={158,66,200},
          thickness=0.5,
          smooth=Smooth.None));
      connect(blood_erythrocytes.solution, others_E.solution) annotation (Line(
          points={{0,-95.4},{0,-88},{72,-88}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(blood_plasma.solution, others_P.solution) annotation (Line(
          points={{0,6.52},{74,6.52},{74,14}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(clock.y, CO2_gas.partialPressure) annotation (Line(
          points={{-33,72},{-24,72},{-24,98},{-84,98},{-84,94}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(H_E.solution, blood_erythrocytes.solution) annotation (Line(
          points={{44.4,-84},{44,-84},{44,-88},{0,-88},{0,-95.4}},
          color={158,66,200},
          smooth=Smooth.None));
        connect(CO2_E.port_a, diffusion.port_b) annotation (Line(
            points={{-70,-72},{-66,-72},{-66,-26}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        connect(CO2.port_a, diffusion.port_a) annotation (Line(
            points={{-68,38},{-66,38},{-66,-6}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        connect(blood_plasma.solution, H.solution) annotation (Line(
            points={{0,6.52},{36,6.52},{36,38},{36.4,38}},
            color={158,66,200},
            smooth=Smooth.None));
        connect(CO2.port_a, HendersonHasselbalch2.substrates[2]) annotation (Line(
            points={{-68,38},{-34,38},{-34,38},{-26,38}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2O.port_a, HendersonHasselbalch2.substrates[1]) annotation (Line(
            points={{-40,22},{-34,22},{-34,34},{-26,34}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        connect(HendersonHasselbalch2.products[1], HCO3.port_a) annotation (Line(
            points={{-6,34},{2,34},{2,24},{8,24}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        connect(HendersonHasselbalch2.products[2], H.port_a) annotation (Line(
            points={{-6,38},{2,38},{2,48},{22,48},{22,47}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.None));
        annotation ( Documentation(info="<html>
<p>CO2 in blood with linear H+ non-bicarbonates buffering without binding to hemoglobin.</p>
<p>The buffer values 0.063 mmol/L commes from Siggaard-Andersen.</p>
</html>",      revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(
          StopTime=1000,
          __Dymola_fixedstepsize=1e-005,
          __Dymola_Algorithm="Lsodar"));
      end CarbonDioxideInBlood;

      model AcidBaseBufferTest
          extends Modelica.Icons.Example;

        Chemical.Sources.Buffer buffer(
          substanceData(z=1.045),
          a_start=10^(-7.2),
          BufferValue=3)
          annotation (Placement(transformation(extent={{-50,4},{-30,24}})));
        Chemical.Components.SimpleSolution simpleSolution
          annotation (Placement(transformation(extent={{-104,-100},{96,100}})));
        Chemical.Sources.ExternalMoleFraction externalMoleFraction(substanceData=
              Chemical.Examples.Substances.Proton_aqueous, MoleFraction=10^(-7.1))
          annotation (Placement(transformation(extent={{0,-46},{20,-26}})));
        Chemical.Components.Substance substance(substanceData=Chemical.Examples.Substances.Water_liquid,
            amountOfSubstance_start=1)
          annotation (Placement(transformation(extent={{52,-82},{72,-62}})));
      equation
        connect(buffer.solution, simpleSolution.solution) annotation (Line(
            points={{-46,4},{-26,4},{-26,-98},{-4,-98}},
            color={158,66,200},
            smooth=Smooth.None));
        connect(externalMoleFraction.port_a, buffer.port_a) annotation (Line(
            points={{20,-36},{40,-36},{40,10},{-30,10},{-30,14}},
            color={158,66,200},
            smooth=Smooth.None));
        connect(simpleSolution.solution, substance.solution) annotation (Line(
            points={{-4,-98},{26,-98},{26,-82},{56,-82}},
            color={158,66,200},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics),
                  experiment(StopTime=0.05));
      end AcidBaseBufferTest;
    end AcidBase;

    package Hemoglobin "Hemoglobin blood gases binding"
      model Allosteric_Hemoglobin_MWC "Monod,Wyman,Changeux (1965)"
        extends Modelica.Icons.Example;

        constant Modelica.SIunits.AmountOfSubstance THb = 0.001
        "Total amount of hemoglobin";

        constant Modelica.SIunits.Temperature T=298.15 "Base Temperature";
        constant Real RT=Modelica.Constants.R*T;

        constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        "Amount of solution used for molarity to mole fraction conversion";

        constant Modelica.SIunits.Volume OneLiter = 0.001;

        constant Real L=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
        constant Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
        constant Modelica.SIunits.Concentration KR=0.000671946
        "Oxygen dissociation coefficient on relaxed(R) hemoglobin subunit";

        constant Real KRx=KR*OneLiter/AmountOfSolutionIn1L
        "Mole fraction based KR";

      //Relative Gibbs formation energies of the substances in the system:
        constant Modelica.SIunits.MolarEnergy
          GO2aq=-RT*log(0.0013/55.508),
          GR0=0,                            GT0=GR0 -RT*log(L),
          GR1=GR0+GO2aq +RT*log(KRx/4),     GT1=GR1 -RT*log(c*L),
          GR2=GR1+GO2aq +RT*log(KRx/(3/2)), GT2=GR2 -RT*log(c^2*L),
          GR3=GR2+GO2aq +RT*log(KRx/(2/3)), GT3=GR3 -RT*log(c^3*L),
          GR4=GR3+GO2aq +RT*log(KRx*4),     GT4=GR4 -RT*log(c^4*L);
                                        //*0.018),

        parameter Real KC = 0.001 "Slow down factor";

      Components.SimpleSolution solution
        annotation (Placement(transformation(extent={{-66,-102},{100,124}})));

        Components.Substance oxygen_unbound(substanceData( DfG_25degC_1bar=GO2aq),
            amountOfSubstance_start(displayUnit="mol") = 1e-5)
          annotation (Placement(transformation(extent={{-62,-46},{-42,-26}})));

        Components.Substance T0(substanceData( DfG_25degC_1bar=GT0), amountOfSubstance_start=
              THb)
          annotation (Placement(transformation(extent={{34,78},{54,98}})));

        Components.Substance T1(substanceData( DfG_25degC_1bar=GT1),
            amountOfSubstance_start=THb*1e-4)
          annotation (Placement(transformation(extent={{34,36},{54,56}})));

        Components.Substance T2(substanceData( DfG_25degC_1bar=GT2),
            amountOfSubstance_start=THb*1e-8)
          annotation (Placement(transformation(extent={{34,-10},{54,10}})));

        Components.Substance R1(substanceData( DfG_25degC_1bar=GR1),
            amountOfSubstance_start=THb*1e-8)
          annotation (Placement(transformation(extent={{-20,36},{0,56}})));

        Components.Substance R2(substanceData( DfG_25degC_1bar=GR2),
            amountOfSubstance_start=THb*1e-10)
          annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

        Components.Substance T3(substanceData( DfG_25degC_1bar=GT3),
            amountOfSubstance_start=THb*1e-12)
          annotation (Placement(transformation(extent={{34,-54},{54,-34}})));

        Components.Substance R3(substanceData( DfG_25degC_1bar=GR3),
            amountOfSubstance_start=THb*1e-12)
          annotation (Placement(transformation(extent={{-20,-54},{0,-34}})));

        Components.Substance T4(substanceData( DfG_25degC_1bar=GT4),
            amountOfSubstance_start=THb*1e-17)
          annotation (Placement(transformation(extent={{34,-92},{54,-72}})));

        Components.Substance R4(substanceData( DfG_25degC_1bar=GR4),
            amountOfSubstance_start=THb*1e-14)
          annotation (Placement(transformation(extent={{-20,-92},{0,-72}})));

        Components.Substance R0(substanceData( DfG_25degC_1bar=GR0),
            amountOfSubstance_start=THb*1e-7)
          annotation (Placement(transformation(extent={{-20,78},{0,98}})));

        Components.Reaction quaternaryForm(KC=KC)
                                           annotation (Placement(transformation(extent={{4,78},{24,98}})));
        Components.Reaction oxyR1(nP=2, KC=KC)
                                        annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-8,64})));
        Components.Reaction oxyT1(nP=2, KC=KC)
                                        annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,64})));
        Components.Reaction oxyR2(nP=2, KC=KC)
                                        annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,22})));
        Components.Reaction oxyR3(nP=2, KC=KC)
                                        annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-24})));
        Components.Reaction oxyR4(nP=2, KC=KC)
                                        annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-10,-66})));
        Components.Reaction oxyT2(nP=2, KC=KC)
                                        annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,22})));
        Components.Reaction oxyT3(nP=2, KC=KC)
                                        annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-24})));
        Components.Reaction oxyT4(nP=2, KC=KC)
                                        annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,-66})));
        Components.Reaction quaternaryForm1(KC=KC)
                                            annotation (Placement(transformation(extent={{8,36},{28,56}})));
        Components.Reaction quaternaryForm2(KC=KC)
                                            annotation (Placement(transformation(extent={{8,-10},{28,10}})));
        Components.Reaction quaternaryForm3(KC=KC)
                                            annotation (Placement(transformation(extent={{8,-54},{28,-34}})));
        Components.Reaction quaternaryForm4(KC=KC)
                                            annotation (Placement(transformation(extent={{10,-92},{30,-72}})));

        Modelica.Blocks.Sources.Clock clock(offset=10)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,62})));
      Sources.ExternalIdealGasSubstance O2_in_air(
        TotalPressure(displayUnit="kPa") = 101325.0144354,
        substanceData=Chemical.Examples.Substances.Oxygen_gas,
        PartialPressure(displayUnit="kPa") = 1000,
        usePartialPressureInput=true) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-84,22})));

        Components.GasSolubility gasSolubility(useWaterCorrection=false, KC=KC)
          annotation (Placement(transformation(extent={{-94,-16},{-74,4}})));

        Components.Substance H2O(substanceData=Chemical.Examples.Substances.Water_liquid,
          amountOfSubstance_start=38.7)
        annotation (Placement(transformation(extent={{-58,-88},{-38,-68}})));

        Real sO2;
      equation
        sO2 = (R1.x + 2*R2.x + 3*R3.x + 4*R4.x + T1.x + 2*T2.x + 3*T3.x + 4*T4.x) /
         (4*(R0.x + R1.x + R2.x + R3.x + R4.x + T0.x + T1.x + T2.x + T3.x + T4.x));

        connect(quaternaryForm.products[1],T0. port_a) annotation (Line(
            points={{24,88},{54,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR1.substrates[1],R1. port_a) annotation (Line(
            points={{-8,54},{-8,46},{0,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R1.port_a,oxyR2. products[1]) annotation (Line(
            points={{0,46},{0,32},{-10.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR2.substrates[1],R2. port_a) annotation (Line(
            points={{-10,12},{-10,0},{0,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.substrates[1],R3. port_a) annotation (Line(
            points={{-10,-34},{-10,-44},{0,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.products[1],R2. port_a) annotation (Line(
            points={{-10.5,-14},{-10.5,-7},{0,-7},{0,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R3.port_a,oxyR4. products[1]) annotation (Line(
            points={{0,-44},{0,-56},{-10.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR4.substrates[1],R4. port_a) annotation (Line(
            points={{-10,-76},{-10,-82},{0,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT1.products[1],T0. port_a) annotation (Line(
            points={{44.5,74},{44.5,88},{54,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT1.substrates[1],T1. port_a) annotation (Line(
            points={{44,54},{44,46},{54,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(T1.port_a,oxyT2. products[1]) annotation (Line(
            points={{54,46},{54,32},{44.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT3.substrates[1],T3. port_a) annotation (Line(
            points={{44,-34},{44,-44},{54,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(T3.port_a,oxyT4. products[1]) annotation (Line(
            points={{54,-44},{54,-56},{44.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT4.substrates[1],T4. port_a) annotation (Line(
            points={{44,-76},{44,-82},{54,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,quaternaryForm. substrates[1]) annotation (Line(
            points={{0,88},{4,88}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0.port_a,oxyR1. products[1]) annotation (Line(
            points={{0,88},{0,74},{-8.5,74}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R1.port_a,quaternaryForm1. substrates[1]) annotation (Line(
            points={{0,46},{8,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm1.products[1],T1. port_a) annotation (Line(
            points={{28,46},{54,46}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R2.port_a,quaternaryForm2. substrates[1]) annotation (Line(
            points={{0,0},{8,0}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R3.port_a,quaternaryForm3. substrates[1]) annotation (Line(
            points={{0,-44},{8,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm3.products[1],T3. port_a) annotation (Line(
            points={{28,-44},{54,-44}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R4.port_a,quaternaryForm4. substrates[1]) annotation (Line(
            points={{0,-82},{10,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm4.products[1],T4. port_a) annotation (Line(
            points={{30,-82},{54,-82}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR1.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-7.5,74},{-42,74},{-42,-36}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR2.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,32},{-42,32},{-42,-36}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR3.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,-14},{-42,-14},{-42,-36}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyR4.products[2],oxygen_unbound. port_a)
                                            annotation (Line(
            points={{-9.5,-56},{-42,-56},{-42,-36}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT1.products[2])
                                            annotation (Line(
            points={{-42,-36},{-42,74},{43.5,74}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT2.products[2])
                                            annotation (Line(
            points={{-42,-36},{-42,32},{43.5,32}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT3.products[2])
                                            annotation (Line(
            points={{-42,-36},{-42,-14},{43.5,-14}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.port_a, oxyT4.products[2])
                                            annotation (Line(
            points={{-42,-36},{-42,-56},{43.5,-56}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(O2_in_air.port_a, gasSolubility.gas_port) annotation (Line(
            points={{-84,12},{-84,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(gasSolubility.liquid_port, oxygen_unbound.port_a) annotation (Line(
            points={{-84,-16},{-84,-36},{-42,-36}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygen_unbound.solution, solution.solution) annotation (Line(
            points={{-58,-46},{-58,-99.74},{17,-99.74}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(R0.solution, solution.solution) annotation (Line(
            points={{-16,78},{-16,-99.74},{17,-99.74}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(T0.solution, solution.solution) annotation (Line(
            points={{38,78},{38,-99.74},{17,-99.74}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(R1.solution, solution.solution) annotation (Line(points={{-16,36},
              {-16,-99.74},{17,-99.74}},
                                  smooth=Smooth.None));
        connect(T1.solution, solution.solution) annotation (Line(points={{38,36},
              {38,-99.74},{17,-99.74}},
                            smooth=Smooth.None));
        connect(R2.solution, solution.solution) annotation (Line(points={{-16,-10},
              {-16,-99.74},{17,-99.74}},
                                  smooth=Smooth.None));
        connect(T3.solution, solution.solution) annotation (Line(points={{38,-54},
              {38,-99.74},{17,-99.74}},
                                  smooth=Smooth.None));
        connect(R3.solution, solution.solution) annotation (Line(points={{-16,-54},
              {-16,-99.74},{17,-99.74}},
                                  smooth=Smooth.None));
        connect(R4.solution, solution.solution) annotation (Line(points={{-16,-92},
              {-16,-99.74},{17,-99.74}},
                                  smooth=Smooth.None));
        connect(T4.solution, solution.solution) annotation (Line(points={{38,-92},
              {38,-99.74},{17,-99.74}},
                                  smooth=Smooth.None));
        connect(quaternaryForm2.products[1], T2.port_a) annotation (Line(
            points={{28,0},{54,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(oxyT2.substrates[1], T2.port_a) annotation (Line(
            points={{44,12},{44,0},{54,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(T2.port_a, oxyT3.products[1]) annotation (Line(
            points={{54,0},{54,-14},{44.5,-14}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(T2.solution, solution.solution) annotation (Line(points={{38,-10},
              {38,-99.74},{17,-99.74}},
                                  smooth=Smooth.None));
        connect(clock.y, O2_in_air.partialPressure) annotation (Line(
            points={{-84,51},{-84,32}},
            color={0,0,127},
            smooth=Smooth.None));
      connect(H2O.solution, solution.solution) annotation (Line(
          points={{-54,-88},{-54,-99.74},{17,-99.74}},
          color={158,66,200},
          smooth=Smooth.None));
        annotation (          experiment(StopTime=15000, Tolerance=0.01),
                                Documentation(info="<html>
<p>To understand the model is necessary to study the principles of MWC allosteric transitions first published by </p>
<p>[1] Monod,Wyman,Changeux (1965). &QUOT;On the nature of allosteric transitions: a plausible model.&QUOT; Journal of molecular biology 12(1): 88-118.</p>
<p><br>In short it is about binding oxygen to hemoglobin.</p>
<p>Oxgen are driven by its partial pressure using clock source - from very little pressure to pressure of 10kPa.</p>
<p>(Partial pressure of oxygen in air is the air pressure multiplied by the fraction of the oxygen in air.)</p>
<p>Hemoglobin was observed (by Perutz) in two structuraly different forms R and T.</p>
<p>These forms are represented by blocks T0..T4 and R0..R4, where the suffexed index means the number of oxygen bounded to the form.</p>
<p><br>In equilibrated model can be four chemical reactions removed and the results will be the same, but dynamics will change a lot. ;)</p>
<p>If you remove the quaternaryForm1,quaternaryForm2,quaternaryForm3,quaternaryForm4 then the model in equilibrium will be exactly the same as in MWC article.</p>
<p><br>Parameters was fitted to data of Severinghaus article from 1979. (For example at pO2=26mmHg is oxygen saturation sO2 = 48.27 &percnt;).</p>
</html>",   revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end Allosteric_Hemoglobin_MWC;

      model Allosteric_Hemoglobin2_MWC
      "Monod,Wyman,Changeux (1965) - The same allosteric hemoglobin model as Allosteric_Hemoglobin_MWC implemented by Speciation blocks"
        extends Modelica.Icons.Example;

        constant Modelica.SIunits.AmountOfSubstance THb = 0.001
        "Total amount of hemoglobin";

        constant Modelica.SIunits.Temperature T=298.15 "Base Temperature";
        constant Real RT=Modelica.Constants.R*T;

        constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        "Amount of solution used for molarity to mole fraction conversion";

        constant Modelica.SIunits.Volume OneLiter = 0.001;

        parameter Real L=7.0529*10^6
        "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
        parameter Real c=0.00431555
        "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
        parameter Modelica.SIunits.Concentration KR=0.000671946
        "oxygen dissociation on relaxed(R) hemoglobin subunit";
                                                                    //*7.875647668393782383419689119171e-5
                                                                  //10.500001495896 7.8756465463794e-05

        parameter Modelica.SIunits.Concentration KT=KR/c
        "oxygen dissociation on tensed(T) hemoglobin subunit";

        parameter Modelica.SIunits.MoleFraction KRx = KR*OneLiter/AmountOfSolutionIn1L;
        parameter Modelica.SIunits.MoleFraction KTx = KT*OneLiter/AmountOfSolutionIn1L;

        parameter Modelica.SIunits.ChemicalPotential DfG_O2 = -RT*log(0.0013/55.508);
        parameter Modelica.SIunits.ChemicalPotential DfG_uR = 0;
        parameter Modelica.SIunits.ChemicalPotential DfG_uRO2 = DfG_uR + DfG_O2 + RT * log(KRx);
        parameter Modelica.SIunits.ChemicalPotential DfG_uT = 0;
        parameter Modelica.SIunits.ChemicalPotential DfG_uTO2 = DfG_uT + DfG_O2 + RT * log(KTx);
        parameter Modelica.SIunits.ChemicalPotential DfG_tT = 0;
        parameter Modelica.SIunits.ChemicalPotential DfG_tR = DfG_tT + RT * log(L);

        parameter Real KC = 1e-3 "Slow down factor";
                                 //0.000001

      Components.SimpleSolution solution
        annotation (Placement(transformation(extent={{-100,-100},{100,58}})));

        Components.Reaction quaternaryForm(KC=KC)
          annotation (Placement(transformation(extent={{-4,-68},{16,-48}})));
        Components.Speciation R0_in_R(NumberOfSubunits=4)
          annotation (Placement(transformation(extent={{-50,-68},{-30,-48}})));
         // AmountOfSubstance_start=4e-11)
        Components.Speciation T0_in_T(NumberOfSubunits=4)
          annotation (Placement(transformation(extent={{68,-68},{48,-48}})));
         // AmountOfSubstance_start=totalAmountOfHemoglobin)
        Components.Substance OxyRHm[4](
              each amountOfSubstance_start=5.88e-9,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG_25degC_1bar=DfG_O2 + RT*log(KRx) + DfG_tR/4))
        "Oxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-96,-18},{-76,2}})));

        Components.Reaction oxygenation_R[4](each nP=2, each KC=KC)
          annotation (Placement(transformation(extent={{-68,-18},{-48,2}})));
        Components.Substance DeoxyRHm[4](
              each amountOfSubstance_start = 1.58e-7,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG_25degC_1bar=DfG_tR/4))
        "Deoxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{0,-34},{-20,-14}})));

        Components.Substance OxyTHm[4](
              each amountOfSubstance_start=1e-4,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG_25degC_1bar=DfG_O2 + RT*log(KTx) + DfG_tT/4))
        "Oxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{14,-18},{34,2}})));

        Components.Reaction oxygenation_T[4](each nP=2, each KC=KC)
          annotation (Placement(transformation(extent={{42,-18},{62,2}})));
        Components.Substance DeoxyTHm[4](
             each  amountOfSubstance_start = THb - 1e-4 - 1.58e-7 - 5.88e-9,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG_25degC_1bar=DfG_tT/4))
        "Deoxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{96,-34},{76,-14}})));

        Components.Substance oxygen_unbound(        substanceData(DfG_25degC_1bar=DfG_O2),
            amountOfSubstance_start=2e-8)
          annotation (Placement(transformation(extent={{-2,6},{18,26}})));
        Modelica.Blocks.Sources.Clock clock(offset=10)
          annotation (Placement(transformation(extent={{-40,74},{-20,94}})));
      Sources.ExternalIdealGasSubstance oxygen_in_air(usePartialPressureInput=
            true, substanceData=Chemical.Examples.Substances.Oxygen_gas)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={8,68})));
        Components.GasSolubility partialPressure1(       useWaterCorrection=false, KC=KC)
                                                  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                origin={8,40})));

        Real sO2 "Hemoglobin oxygen saturation";
        Components.Substance H2O(substanceData=Chemical.Examples.Substances.Water_liquid,
            amountOfSubstance_start=38.7)
          annotation (Placement(transformation(extent={{68,-94},{88,-74}})));
      equation
        sO2 = (sum(OxyRHm.x) + sum(OxyTHm.x)) /
        (sum(DeoxyRHm.x) + sum(DeoxyTHm.x) + sum(OxyRHm.x) + sum(OxyTHm.x));

        connect(OxyTHm.port_a, oxygenation_T.substrates[1])
                                                 annotation (Line(
            points={{34,-8},{42,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_T.products[1], DeoxyTHm.port_a)
                                               annotation (Line(
            points={{62,-8.5},{66,-8.5},{66,-24},{76,-24}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(clock.y, oxygen_in_air.partialPressure) annotation (Line(
            points={{-19,84},{8,84},{8,78}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(OxyRHm.port_a, oxygenation_R.substrates[1]) annotation (Line(
            points={{-76,-8},{-68,-8}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(DeoxyRHm.port_a, R0_in_R.subunits) annotation (Line(
            points={{-20,-24},{-42,-24},{-42,-48}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(oxygenation_R.products[1], DeoxyRHm.port_a) annotation (Line(
            points={{-48,-8.5},{-34,-8.5},{-34,-24},{-20,-24}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(T0_in_T.subunits, DeoxyTHm.port_a)   annotation (Line(
            points={{60,-48},{60,-24},{76,-24}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));

        connect(oxygen_in_air.port_a, partialPressure1.gas_port) annotation (Line(
            points={{8,58},{8,50}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(partialPressure1.liquid_port, oxygen_unbound.port_a) annotation (Line(
            points={{8,30},{8,16},{18,16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(R0_in_R.port_a, quaternaryForm.substrates[1]) annotation (Line(
            points={{-30,-68},{-18,-68},{-18,-58},{-4,-58}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(quaternaryForm.products[1], T0_in_T.port_a) annotation (Line(
            points={{16,-58},{32,-58},{32,-68},{48,-68}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));

        for i in 1:4 loop
          connect(oxygenation_T[i].products[2], oxygen_unbound.port_a) annotation (Line(
            points={{62,-7.5},{70,-7.5},{70,16},{18,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
          connect(oxygenation_R[i].products[2], oxygen_unbound.port_a) annotation (Line(
            points={{-48,-7.5},{-12,-7.5},{-12,16},{18,16}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.None));
        connect(R0_in_R.subunitSolution, DeoxyRHm[i].solution) annotation (Line(
            points={{-36,-48},{-36,-42},{-4,-42},{-4,-34}},
            color={158,66,200},
            smooth=Smooth.None));
        connect(R0_in_R.subunitSolution, OxyRHm[i].solution) annotation (Line(
            points={{-36,-48},{-36,-42},{-92,-42},{-92,-18}},
            color={158,66,200},
            smooth=Smooth.None));
        connect(OxyTHm[i].solution, T0_in_T.subunitSolution) annotation (Line(
            points={{18,-18},{18,-44},{54,-44},{54,-48}},
            color={158,66,200},
            smooth=Smooth.None));
        connect(DeoxyTHm[i].solution, T0_in_T.subunitSolution) annotation (Line(
            points={{92,-34},{92,-44},{54,-44},{54,-48}},
            color={158,66,200},
            smooth=Smooth.None));
        end for;

        connect(R0_in_R.solution, solution.solution) annotation (Line(
            points={{-46,-68},{-46,-98.42},{0,-98.42}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(T0_in_T.solution, solution.solution) annotation (Line(
            points={{64,-68},{64,-98.42},{0,-98.42}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(oxygen_unbound.solution, solution.solution) annotation (Line(points={{2,6},{2,
                -98.42},{0,-98.42}},             smooth=Smooth.None));
        connect(solution.solution, H2O.solution) annotation (Line(
            points={{0,-98.42},{72,-98.42},{72,-94}},
            color={158,66,200},
            smooth=Smooth.None));

        annotation (          experiment(StopTime=15000, __Dymola_Algorithm="Dassl"),
          Documentation(revisions="<html>
<p><i>2013-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p>Before silumation in &QUOT;Dymola 2014 FD01&QUOT; please chose &QUOT;Euler&QUOT; method!</p>
<p><br>To understand the model is necessary to study the principles of MWC allosteric transitions first published by </p>
<p>[1] Monod,Wyman,Changeux (1965). &QUOT;On the nature of allosteric transitions: a plausible model.&QUOT; Journal of molecular biology 12(1): 88-118.</p>
<p><br>In short it is about binding oxygen to hemoglobin.</p>
<p>Oxgen are driven by its partial pressure using clock source - from very little pressure to pressure of 10kPa.</p>
<p>(Partial pressure of oxygen in air is the air pressure multiplied by the fraction of the oxygen in air.)</p>
<p>Hemoglobin was observed (by Perutz) in two structuraly different forms R and T.</p>
<p>These forms are represented by blocks T0..T4 and R0..R4, where the suffexed index means the number of oxygen bounded to the form.</p>
<p><br>In equilibrated model can be four chemical reactions removed and the results will be the same, but dynamics will change a lot. ;)</p>
<p>If you remove the quaternaryForm1,quaternaryForm2,quaternaryForm3,quaternaryForm4 then the model in equilibrium will be exactly the same as in MWC article.</p>
<p><br>Parameters was fitted to data of Severinghaus article from 1979. (For example at pO2=26mmHg is oxygen saturation sO2 = 48.27 &percnt;).</p>
</html>"));
      end Allosteric_Hemoglobin2_MWC;
    end Hemoglobin;

    package CheckSubstancesData
      model SimpleReaction
      "The simple chemical reaction A<->B with equilibrium B/A = 2"
         extends Modelica.Icons.Example;

        constant Real K = 2 "Dissociation constant of the reaction";

        constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
        constant Real R = Modelica.Constants.R "Gas constant";

      Sensors.DissociationCoefficient dissociationCoefficient
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Sources.PureSubstance A
        annotation (Placement(transformation(extent={{-56,-10},{-36,10}})));
      Sources.PureSubstance B(redeclare package stateOfMatter =
            Chemical.Interfaces.Incompressible, substanceData(DfG_25degC_1bar=-
              R*T_25degC*log(K)))
        annotation (Placement(transformation(extent={{60,-10},{40,10}})));
      equation
      connect(A.port_a, dissociationCoefficient.substrates[1]) annotation (Line(
          points={{-36,0},{-10,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
      connect(B.port_a, dissociationCoefficient.products[1]) annotation (Line(
          points={{40,0},{10,0}},
          color={158,66,200},
          thickness=1,
          smooth=Smooth.None));
        annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
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

        constant Modelica.SIunits.Temperature T_25degC = 298.15 "Temperature";
        constant Real R = Modelica.Constants.R "Gas constant";

      Sources.PureSubstance A
        annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
        Sensors.DissociationCoefficient
                            reaction(nS=2)
          annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Sources.PureSubstance B
        annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Sources.PureSubstance C(substanceData(DfG_25degC_1bar=-R*T_25degC*log(Kx)))
        annotation (Placement(transformation(extent={{68,-8},{48,12}})));

      equation
        connect(reaction.products[1], C.port_a) annotation (Line(
            points={{24,2},{48,2}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));

        connect(B.port_a, reaction.substrates[1]) annotation (Line(
            points={{-14,-14},{-10,-14},{-10,1.5},{4,1.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(A.port_a, reaction.substrates[2]) annotation (Line(
            points={{-14,12},{-10,12},{-10,2.5},{4,2.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.001));
      end SimpleReaction2;

      model SimpleReaction2_Get_DfG
      "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
         extends Modelica.Icons.Example;

      Sources.PureSubstance A
        annotation (Placement(transformation(extent={{-28,42},{-8,62}})));
        Sensors.DissociationCoefficient
                            reaction(nS=2)
          annotation (Placement(transformation(extent={{10,32},{30,52}})));
      Sources.PureSubstance B
        annotation (Placement(transformation(extent={{-28,16},{-8,36}})));

        Modelica.Blocks.Math.InverseBlockConstraints inverseBlockConstraints
          annotation (Placement(transformation(extent={{-42,-80},{82,80}})));
      Sources.ExternalElectroChemicalPotential C(usePotentialInput=true)
        annotation (Placement(transformation(extent={{60,32},{40,52}})));
        Modelica.Blocks.Sources.Constant K(k=2*55.508)
          annotation (Placement(transformation(extent={{-92,-10},{-72,10}})));
      equation

        connect(B.port_a, reaction.substrates[1]) annotation (Line(
            points={{-8,26},{-4,26},{-4,41.5},{10,41.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(A.port_a, reaction.substrates[2]) annotation (Line(
            points={{-8,52},{-4,52},{-4,42.5},{10,42.5}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(reaction.DissociationCoefficient_MoleFractionBased,
          inverseBlockConstraints.u2) annotation (Line(
            points={{20,34},{20,0},{-29.6,0}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(reaction.products[1], C.port_a) annotation (Line(
            points={{30,42},{40,42}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(C.uInput, inverseBlockConstraints.y2) annotation (Line(
            points={{60,42},{70,42},{70,24},{46,24},{46,0},{72.7,0}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(inverseBlockConstraints.u1, K.y) annotation (Line(
            points={{-48.2,0},{-71,0}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation ( Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.001));
      end SimpleReaction2_Get_DfG;

      model StandardElectrochemicalCell
      "Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential "
       extends Modelica.Icons.Example;

      Components.SimpleSolution cathode(ElectricGround=false)
        annotation (Placement(transformation(extent={{-90,-40},{-46,68}})));

      Components.SimpleSolution anode(ElectricGround=false)
        annotation (Placement(transformation(extent={{60,-40},{96,70}})));

      Sources.PureSubstance Ag(substanceData=Chemical.Examples.Substances.Silver_solid)
        annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
      Sources.PureSubstance Cl(substanceData=Chemical.Examples.Substances.Chloride_aqueous)
        annotation (Placement(transformation(extent={{-8,-36},{-28,-16}})));
      Sources.PureSubstance AgCl(substanceData=Chemical.Examples.Substances.SilverChloride_solid)
        annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
      Sources.ExternalIdealGasSubstance H2(
        substanceData=Chemical.Examples.Substances.Hydrogen_gas,
        PartialPressure=100000,
        TotalPressure=100000)
        annotation (Placement(transformation(extent={{24,32},{44,52}})));
      Sources.PureSubstance H(substanceData=Chemical.Examples.Substances.Proton_aqueous)
        annotation (Placement(transformation(extent={{18,-36},{38,-16}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-6,64},{14,84}})));
        Components.Reaction electrodeReaction(nP=2, p={2,2}) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={52,6})));
        Components.Reaction electrodeReaction1(nS=2, nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-40,6})));
      Components.ElectronTransfer electrone
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
      Components.ElectronTransfer electrone1
        annotation (Placement(transformation(extent={{86,-26},{66,-6}})));

      equation
        connect(Ag.port_a, electrodeReaction1.substrates[1]) annotation (Line(
            points={{-60,-18},{-42,-18},{-42,-4},{-40.5,-4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(Cl.port_a, electrodeReaction1.substrates[2]) annotation (Line(
            points={{-28,-26},{-39.5,-26},{-39.5,-4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(AgCl.port_a, electrodeReaction1.products[1]) annotation (Line(
            points={{-60,22},{-40.5,22},{-40.5,16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H2.port_a, electrodeReaction.substrates[1]) annotation (Line(
            points={{44,42},{52,42},{52,16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H.port_a, electrodeReaction.products[1]) annotation (Line(
            points={{38,-26},{52.5,-26},{52.5,-4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(electrodeReaction.products[2], electrone1.port_a) annotation (Line(
            points={{51.5,-4},{51.5,-16},{66,-16}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(electrodeReaction1.products[2], electrone.port_a) annotation (Line(
            points={{-39.5,16},{-39.5,50},{-60,50}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(electrone.pin, voltageSensor.p) annotation (Line(
            points={{-80,50},{-86,50},{-86,74},{-6,74}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(electrone1.pin, voltageSensor.n) annotation (Line(
            points={{86,-16},{90,-16},{90,74},{14,74}},
            color={0,0,255},
            smooth=Smooth.None));
      connect(electrone1.solution, anode.solution) annotation (Line(
          points={{82,-26},{80,-26},{80,-38.9},{78,-38.9}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(electrone.solution, cathode.solution) annotation (Line(
          points={{-76,40},{-88,40},{-88,-38.92},{-68,-38.92}},
          color={158,66,200},
          smooth=Smooth.None));
        annotation (
        experiment(StopTime=1), Documentation(info=
                      "<html>
<p>Hypothetical experiment of pure substances reaction to define the standard electrochemical cell potential </p>
</html>",   revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end StandardElectrochemicalCell;

      model StandardLeadAcidPotential
      "Standard potential of the lead acid battery"
       extends Modelica.Icons.Example;

      Components.SimpleSolution anode(ElectricGround=
            false)
        annotation (Placement(transformation(extent={{54,-46},{92,62}})));

      Components.SimpleSolution cathode(ElectricGround=false)
        annotation (Placement(transformation(extent={{-94,-50},{-56,58}})));

      Sources.PureSubstance Pb(substanceData=Chemical.Examples.Substances.Lead_solid)
        annotation (Placement(transformation(extent={{84,-34},{64,-14}})));
      Sources.PureSubstance HSO4(substanceData=Chemical.Examples.Substances.HydrogenSulfate_aqueous)
        annotation (Placement(transformation(extent={{-22,-58},{-2,-38}})));
      Sources.PureSubstance PbSO4_(substanceData=Chemical.Examples.Substances.LeadSulfate_solid)
        annotation (Placement(transformation(extent={{84,4},{64,24}})));
      Sources.PureSubstance H(substanceData=Chemical.Examples.Substances.Proton_aqueous)
        annotation (Placement(transformation(extent={{6,-28},{26,-8}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-6,60},{14,80}})));
        Components.Reaction electrodeReaction(nP=2,
          nS=4,
          s={1,1,3,2},
        p={1,2})                                             annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-42,14})));
        Components.Reaction electrodeReaction1(nS=2,
          nP=3,
          p={1,1,2})                                       annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={44,14})));

      Components.ElectronTransfer electrone
        annotation (Placement(transformation(extent={{84,32},{64,52}})));
      Components.ElectronTransfer electrone1
        annotation (Placement(transformation(extent={{-86,-12},{-66,8}})));
      Sources.PureSubstance PbO2(substanceData=Chemical.Examples.Substances.LeadDioxide_solid)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            origin={-74,-30})));
      Sources.PureSubstance H2O(substanceData=Chemical.Examples.Substances.Water_liquid)
        annotation (Placement(transformation(extent={{-2,-10},{-22,10}})));
      Sources.PureSubstance PbSO4(substanceData=Chemical.Examples.Substances.LeadSulfate_solid)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            origin={-74,32})));

      equation
        connect(Pb.port_a, electrodeReaction1.substrates[1]) annotation (Line(
            points={{64,-24},{45.5,-24},{45.5,4},{44.5,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HSO4.port_a, electrodeReaction1.substrates[2]) annotation (Line(
            points={{-2,-48},{44,-48},{44,4},{43.5,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(PbSO4_.port_a, electrodeReaction1.products[1]) annotation (Line(
            points={{64,14},{56,14},{56,28},{46,28},{46,24},{44.6667,24}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(electrodeReaction.products[1], PbSO4.port_a) annotation (Line(
            points={{-42.5,24},{-42.5,32},{-64,32}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(electrodeReaction.products[2], H2O.port_a) annotation (Line(
            points={{-41.5,24},{-40,24},{-40,32},{-34,32},{-34,0},{-22,0}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(PbO2.port_a, electrodeReaction.substrates[1]) annotation (Line(
            points={{-64,-30},{-42,-30},{-42,4},{-42.75,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(HSO4.port_a, electrodeReaction.substrates[2]) annotation (Line(
            points={{-2,-48},{-40,-48},{-40,4},{-42.25,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H.port_a, electrodeReaction.substrates[3]) annotation (Line(
            points={{26,-18},{-38,-18},{-38,4},{-41.75,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(electrone1.port_a, electrodeReaction.substrates[4]) annotation (Line(
            points={{-66,-2},{-44,-2},{-44,4},{-41.25,4}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(
            points={{26,-18},{32,-18},{32,32},{42,32},{42,24},{44,24}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
        connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
            points={{64,42},{44,42},{44,24},{43.3333,24}},
            color={158,66,200},
            thickness=1,
            smooth=Smooth.None));
      connect(electrone1.pin, voltageSensor.p) annotation (Line(
          points={{-86,-2},{-98,-2},{-98,70},{-6,70}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(electrone.pin, voltageSensor.n) annotation (Line(
          points={{84,42},{96,42},{96,70},{14,70}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(electrone1.solution, cathode.solution) annotation (Line(
          points={{-82,-12},{-82,-48.92},{-75,-48.92}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(electrone.solution, anode.solution) annotation (Line(
          points={{80,32},{80,-44.92},{73,-44.92}},
          color={158,66,200},
          smooth=Smooth.None));
        annotation (
        experiment(StopTime=100), Documentation(revisions=
                        "<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end StandardLeadAcidPotential;

    end CheckSubstancesData;

    model FluidAdapter
     extends Modelica.Icons.Example;

      Chemical.Components.FluidAdapter fluidConversion(
        redeclare package Medium =
            Modelica.Media.CompressibleLiquids.LinearWater_pT_Ambient,
        substanceNames={"H2O(l)"},
        substanceData={Chemical.Examples.Substances.Water_liquid})
        annotation (Placement(transformation(extent={{-16,26},{4,46}})));
      Chemical.Components.SimpleSolution simpleSolution
        annotation (Placement(transformation(extent={{-74,8},{16,70}})));
      Modelica.Fluid.Pipes.StaticPipe pipe(
        redeclare package Medium =
            Modelica.Media.CompressibleLiquids.LinearWater_pT_Ambient,
        length=1,
        diameter=1) annotation (Placement(transformation(extent={{24,26},{44,46}})));
      Modelica.Fluid.Vessels.OpenTank tank(
        nPorts=1,
        height=1,
        crossArea=1,
        redeclare package Medium =
            Modelica.Media.CompressibleLiquids.LinearWater_pT_Ambient,
        portsData={Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter=1)})
        annotation (Placement(transformation(extent={{52,42},{92,82}})));
      inner Modelica.Fluid.System system
        annotation (Placement(transformation(extent={{-92,76},{-72,96}})));
      Chemical.Components.Substance H2O(
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=55)
        annotation (Placement(transformation(extent={{-58,26},{-38,46}})));
      Chemical.Components.FluidAdapter fluidConversion1(
        redeclare package Medium =
            Modelica.Media.CompressibleLiquids.LinearWater_pT_Ambient,
        substanceNames={"H2O(l)"},
        substanceData={Chemical.Examples.Substances.Water_liquid})
        annotation (Placement(transformation(extent={{-52,-56},{-32,-36}})));
      Chemical.Components.SimpleSolution simpleSolution1
        annotation (Placement(transformation(extent={{-96,-74},{-26,-14}})));
      Chemical.Components.Substance H2O1(
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=55)
        annotation (Placement(transformation(extent={{-80,-56},{-60,-36}})));
      Components.Solution                simpleSolution2(T(start=298))
        annotation (Placement(transformation(extent={{24,-74},{98,-12}})));
      Chemical.Components.Substance H2O2(
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        substanceData=Chemical.Examples.Substances.Water_liquid,
        amountOfSubstance_start=55)
        annotation (Placement(transformation(extent={{84,-56},{64,-36}})));
      Chemical.Components.FluidAdapter fluidConversion2(
        redeclare package Medium =
            Modelica.Media.CompressibleLiquids.LinearWater_pT_Ambient,
        substanceNames={"H2O(l)"},
        substanceData={Chemical.Examples.Substances.Water_liquid})
        annotation (Placement(transformation(extent={{56,-56},{36,-36}})));
      Modelica.Fluid.Pipes.StaticPipe pipe1(
        redeclare package Medium =
            Modelica.Media.CompressibleLiquids.LinearWater_pT_Ambient,
        length=1,
        diameter=1) annotation (Placement(transformation(extent={{-10,-56},{10,
              -36}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed1
        annotation (Placement(transformation(extent={{50,-90},{70,-70}})));
      Modelica.Mechanics.Translational.Components.Spring spring(c=1)   annotation (
          Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={60,2})));
      Modelica.Mechanics.Translational.Components.Fixed fixed(s0=1)
        annotation (Placement(transformation(extent={{76,0},{96,20}})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{4,-32},{24,-12}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature coolerTemperature(T=298.15)
        annotation (Placement(transformation(extent={{-28,-100},{-8,-80}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=2)
        annotation (Placement(transformation(extent={{26,-100},{6,-80}})));
    equation
      connect(fluidConversion.solution, simpleSolution.solution) annotation (Line(
          points={{-6,31},{-6,8.62},{-2,8.62}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(fluidConversion.fluid, pipe.port_a) annotation (Line(
          points={{4,36},{24,36}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(pipe.port_b, tank.ports[1]) annotation (Line(
          points={{44,36},{72,36},{72,42}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(H2O.port_a, fluidConversion.substances[1]) annotation (Line(
          points={{-38,36},{-16,36}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(H2O.solution, simpleSolution.solution) annotation (Line(
          points={{-54,26},{-54,8.62},{-2,8.62}},
          color={0,128,255},
          smooth=Smooth.None));
    connect(fluidConversion1.solution, simpleSolution1.solution) annotation (
        Line(
        points={{-42,-51},{-42,-73.4},{-40,-73.4}},
        color={0,128,255},
        smooth=Smooth.None));
    connect(H2O1.port_a, fluidConversion1.substances[1]) annotation (Line(
        points={{-60,-46},{-52,-46}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(H2O1.solution, simpleSolution1.solution) annotation (Line(
        points={{-76,-56},{-76,-73.4},{-40,-73.4}},
        color={0,128,255},
        smooth=Smooth.None));
    connect(fluidConversion2.solution, simpleSolution2.solution) annotation (
        Line(
        points={{46,-51},{46,-73.38},{83.2,-73.38}},
        color={0,128,255},
        smooth=Smooth.None));
    connect(H2O2.solution, simpleSolution2.solution) annotation (Line(
        points={{80,-56},{80,-73.38},{83.2,-73.38}},
        color={0,128,255},
        smooth=Smooth.None));
    connect(H2O2.port_a, fluidConversion2.substances[1]) annotation (Line(
        points={{64,-46},{56,-46}},
        color={158,66,200},
        smooth=Smooth.None));
    connect(fluidConversion1.fluid, pipe1.port_a) annotation (Line(
        points={{-32,-46},{-10,-46}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(fluidConversion2.fluid, pipe1.port_b) annotation (Line(
        points={{36,-46},{10,-46}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(spring.flange_a, fixed.flange) annotation (Line(
        points={{60,12},{74,12},{74,10},{86,10}},
        color={0,127,0},
        smooth=Smooth.None));
      connect(thermalConductor.port_b,coolerTemperature. port) annotation (Line(
          points={{6,-90},{-8,-90}},
          color={191,0,0},
          smooth=Smooth.None));
    connect(simpleSolution2.surfaceFlange, spring.flange_b) annotation (Line(
        points={{61,-12},{60,-12},{60,-8}},
        color={0,127,0},
        smooth=Smooth.None));
    connect(simpleSolution2.bottom, fixed1.flange) annotation (Line(
        points={{61,-74.62},{60,-74.62},{60,-80}},
        color={0,127,0},
        smooth=Smooth.None));
    connect(ground.p, simpleSolution2.electricPin) annotation (Line(
        points={{14,-12},{38.8,-12}},
        color={0,0,255},
        smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),   graphics));
    end FluidAdapter;
  end Examples;


  package Components "Chemical Components"
    extends Modelica.Icons.Package;
    model SimpleSolution
    "Chemical solution as homogenous mixture of the substances at constant pressure, grounded or electroneutral, isothermal or at constant temperature"
      extends Icons.Solution;

      extends Interfaces.PartialSolution(T(start=temperature_start),p(start=ConstantPressure));

      parameter Modelica.SIunits.Temperature temperature_start=298.15
      "Initial temperature of the solution"
         annotation (Dialog(group="Initialization"));

      parameter Modelica.SIunits.Pressure ConstantPressure=100000
      "Constant pressure of the solution";

      parameter Boolean ElectricGround = true
      "Is the solution electric potential equal to zero during simulation?";

      parameter Boolean ConstantTemperature = true
      "Has the solution constant temperature during simulation?";

    initial equation
      T=temperature_start;
    equation
      //hydraulic
      solution.p = ConstantPressure;
      workFromEnvironment = 0;

      //electric
      if  ElectricGround then
        //Solution connected to ground has zero voltage. However, electric current from the solution can varies.
        solution.v = 0;
      else
        //Electrically isolated solution has not any electric current from/to the solution. However, electric potential can varies.
        solution.i = 0;
      end if;

      //thermal
      if ConstantTemperature then
        //Ideal thermal exchange between environment and solution to reach constant temperature
        der(T)=0;
      else
        //Thermally isolated without any thermal exchange with environment
        heatFromEnvironment = 0;
      end if;

                                                                                                        annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false, initialScale=1, extent={{-100,-100},{
              100,100}}),
            graphics={Text(
              extent={{-92,-86},{76,-94}},
              lineColor={0,0,255},
              textString="%name",
              horizontalAlignment=TextAlignment.Left)}),
        Documentation(revisions="<html>
<p>2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<h4>amountOfSolution = &int; molarFlows</h4>
<h4>mass = &int; massChanges</h4>
<h4>volume = &int; volumeChanges</h4>
<h4>freeEnthalpy = &int; enthalpyChanges</h4>
<h4>freeEntropy = &int; entropyChanges</h4>
<h4>freeGibbsEnergy = &int; freeGibbsEnergyChanges</h4>
<p>Integration of all substances together into one homogenous mixture - the solution.</p>
</html>"));
    end SimpleSolution;

    model Solution "Chemical solution as homogenous mixture of the substances"
      extends Icons.Solution;

      extends Interfaces.PartialSolution(T(start=temperature_start),p(start=AmbientPressure));

      parameter Modelica.SIunits.Area SurfaceArea=0.01
      "Area for surfacePort to connect MultiBody components";

      parameter Modelica.SIunits.Pressure AmbientPressure=100000
      "Ambient pressure if the force on port surfaceFlange is zero";

       parameter Modelica.SIunits.Temperature temperature_start=298.15
      "Initial temperature of the solution"
         annotation (Dialog(group="Initialization"));

       parameter Boolean isPistonPositionAbsolute=false
      "Relavite position has zero at initial state without force";

      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation (
          Placement(transformation(extent={{-70,-90},{-50,-70}}),iconTransformation(
              extent={{-62,-104},{-58,-100}})));
      Modelica.Electrical.Analog.Interfaces.PositivePin electricPin annotation (Placement(
            transformation(extent={{-70,70},{-50,90}}),    iconTransformation(
              extent={{-62,98},{-58,102}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_a surfaceFlange
      "The pressure of solution generate force on prescribed surface."
        annotation (Placement(transformation(extent={{-10,70},{10,90}}),
            iconTransformation(extent={{-2,98},{2,102}})));

  protected
      parameter Modelica.SIunits.Position positionShift(fixed=false)
      "=0 absolute, otherwise negative";
  public
      Modelica.Mechanics.Translational.Interfaces.Flange_b bottom
      "Fix of the cilinder on bottom."   annotation (Placement(transformation(
              extent={{-10,-90},{10,-70}}), iconTransformation(extent={{-2,-104},{2,
                -100}})));
    initial equation
      T=temperature_start;
      positionShift= if
                       (isPistonPositionAbsolute) then 0 else volume/SurfaceArea;
    equation
      //hydraulic
      surfaceFlange.f+bottom.f = 0;
      surfaceFlange.s-bottom.s = volume/SurfaceArea;
      workFromEnvironment = der(surfaceFlange.f*surfaceFlange.s); //=der( (p-p0) * volume)
      solution.p = AmbientPressure + surfaceFlange.f/SurfaceArea;

      //electric
      electricPin.v=solution.v;
      electricPin.i=solution.i;

      //thermal
      heatPort.T = T;
      heatPort.Q_flow = heatFromEnvironment;

                                                                                                        annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false, initialScale=1, extent={{-100,-100},{100,100}}),
            graphics={Text(
              extent={{-92,-86},{76,-94}},
              lineColor={0,0,255},
              textString="%name",
              horizontalAlignment=TextAlignment.Left)}),
        Documentation(revisions="<html>
<p>2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<h4>amountOfSolution = &sum; amountOfSubstances</h4>
<h4>mass = &sum; massOfSubstances</h4>
<h4>volume = &sum; volumeOfSubstances</h4>
<h4>freeGibbsEnergy = &sum; freeGibbsEnergiesOfSubstances</h4>
<p>To calculate the sum of extensive substance's properties is misused the Modelica \"flow\" prefix even there are not real physical flows. </p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end Solution;

    model IdealGasSolution
    "Chemical solution as gaseous homogenous mixture of the substances"
      extends Icons.Solution;

      extends Interfaces.PartialSolution(T(start=temperature_start),p(start=AmbientPressure));
                                        /*(final amountOfSolution_start=AmbientPressure*volume_start/(Modelica.Constants.R*temperature_start),
  final mass_start=MolarMass_start*AmbientPressure*volume_start/(Modelica.Constants.R*temperature_start))*/

      parameter Modelica.SIunits.Area SurfaceArea=0.01
      "Area for surfacePort to connect MultiBody components";

      parameter Modelica.SIunits.Pressure AmbientPressure=100000
      "Ambient pressure if the force on port surfaceFlange is zero";

       parameter Boolean ElectricGround = true
      "Is the solution electric potential equal to zero during simulation?";

      parameter Modelica.SIunits.Temperature temperature_start=298.15
      "Initial temperature of the solution"
         annotation (Dialog(group="Initialization"));

      parameter Boolean isPistonPositionAbsolute=false
      "Relavite position has zero at initial state without force";

      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation (
          Placement(transformation(extent={{-70,-92},{-50,-72}}),iconTransformation(
              extent={{-62,-104},{-58,-100}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_a surfaceFlange
      "The pressure of solution generate force on prescribed surface."
        annotation (Placement(transformation(extent={{-10,70},{10,90}}),
            iconTransformation(extent={{-2,98},{2,102}})));
  protected
      parameter Modelica.SIunits.Position positionShift(fixed=false)
      "=0 absolute, otherwise negative";
  public
      Modelica.Mechanics.Translational.Interfaces.Flange_b bottom
      "Fix of the cilinder on bottom."   annotation (Placement(transformation(
              extent={{-10,-90},{10,-70}}), iconTransformation(extent={{-2,-104},{2,
                -100}})));
    initial equation
      T=temperature_start;
      positionShift= if
                       (isPistonPositionAbsolute) then 0 else volume/SurfaceArea;
    equation
      //hydraulic
      solution.p = AmbientPressure - surfaceFlange.f/SurfaceArea;

      workFromEnvironment = -der(surfaceFlange.f*surfaceFlange.s); //=der( (p-p0) * volume)

      surfaceFlange.f+bottom.f = 0;
      surfaceFlange.s-bottom.s = volume/SurfaceArea - positionShift;

      //electric
       if  ElectricGround then
        //Solution connected to ground has zero voltage. However, electric current from the solution can varies.
        solution.v = 0;
      else
        //Electrically isolated solution has not any electric current from/to the solution. However, electric potential can varies.
        solution.i = 0;
      end if;

      //thermal
      heatPort.T = T;
      heatPort.Q_flow = heatFromEnvironment;

                                                                                                        annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false, initialScale=1, extent={{-100,-100},{100,100}}),
            graphics={Text(
              extent={{-92,-86},{76,-94}},
              lineColor={0,0,255},
              textString="%name",
              horizontalAlignment=TextAlignment.Left)}),
        Documentation(revisions="<html>
<p>2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<h4>amountOfSolution = &int; molarFlows</h4>
<h4>mass = &int; massChanges</h4>
<h4>volume = &int; volumeChanges</h4>
<h4>freeEnthalpy = &int; enthalpyChanges</h4>
<h4>freeEntropy = &int; entropyChanges</h4>
<h4>freeGibbsEnergy = &int; freeGibbsEnergyChanges</h4>
<p>Integration of all substances together into one homogenous mixture - the solution.</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end IdealGasSolution;

    model Substance "Substance in solution"
      extends Icons.Substance;

      Modelica.SIunits.Concentration c "Molar concentration";

      extends Interfaces.PartialSubstanceInSolution;

      //If it is selected the amount of solution per one kilogram of solvent then the values of amountOfSubstance will be the same as molality
      //If it is selected the amount of solution in one liter of solution then the values of amountOfSubstance will be the same as molarity
      parameter Modelica.SIunits.AmountOfSubstance amountOfSubstance_start=1e-8
      "Initial amount of the substance in compartment"   annotation(HideResult=true);

  protected
      Modelica.SIunits.AmountOfSubstance amountOfSubstance(start=amountOfSubstance_start);
      Real log10n(stateSelect=StateSelect.prefer, start=log10(amountOfSubstance_start))
      "Decadic logarithm of the amount of the substance in solution";
      constant Real InvLog_10=1/log(10);

    initial equation

      amountOfSubstance=amountOfSubstance_start;

    equation

      //The main accumulation equation is "der(amountOfSubstance)=port_a.q"
      // However, the numerical solvers can handle it in form of log10n much better. :-)
      der(log10n)=(InvLog_10)*(port_a.q/amountOfSubstance);
      amountOfSubstance = 10^log10n;

      //Molar Concentration
      c = amountOfSubstance/solution.V;

      //Mole fraction is an analogy of molar concentration or molality.
      x = amountOfSubstance/solution.n;

      //solution flows
      solution.dH = molarEnthalpy*port_a.q + der(molarEnthalpy)*amountOfSubstance;
      solution.i = Modelica.Constants.F * z * port_a.q + Modelica.Constants.F*der(z)*amountOfSubstance;
      solution.dV = molarVolume * port_a.q + der(molarVolume)*amountOfSubstance;

      //extensive properties
      solution.nj=amountOfSubstance;
      solution.mj=amountOfSubstance*molarMass;
      solution.Vj=amountOfSubstance*molarVolume;
      solution.Gj=amountOfSubstance*port_a.u;
      solution.Qj=Modelica.Constants.F*amountOfSubstance*z;
      solution.Ij=(1/2) * ( amountOfSubstance * z^2);

                                                                                                        annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={Text(
              extent={{-84,20},{94,60}},
              lineColor={0,0,255},
              textString="%name")}),
        Documentation(revisions="<html>
<p>2009-2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<h4>n = x &middot; n(solution) = &int; MolarFlow</h4>
<p>where n is amount of the substance and x is mole fraction.</p>
<p>The main class from &ldquo;Chemical&rdquo; package is called &QUOT;Substance&QUOT;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &QUOT;n&QUOT; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &QUOT;n(solution)&QUOT; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solution. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
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

    model Reaction "Chemical Reaction"

      parameter Integer nS=1 "Number of substrates types"
        annotation ( HideResult=true);

      parameter Modelica.SIunits.StoichiometricNumber s[nS]=ones(nS)
      "Stoichiometric reaction coefficient for substrates"
        annotation (HideResult=true);

      parameter Integer nP=1 "Number of products types"
        annotation ( HideResult=true);

      parameter Modelica.SIunits.StoichiometricNumber p[nP]=ones(nP)
      "Stoichiometric reaction coefficients for products"
        annotation (HideResult=true);

      Modelica.SIunits.MolarFlowRate rr(start=0) "Reaction molar flow rate";



      extends Interfaces.ConditionalKinetics;

    Interfaces.SubstancePorts_b substrates[nS] annotation (Placement(
          transformation(extent={{-110,-40},{-90,40}}), iconTransformation(
            extent={{-110,-40},{-90,40}})));
    Interfaces.SubstancePorts_b products[nP] annotation (Placement(
          transformation(extent={{90,-40},{110,40}}), iconTransformation(extent
            ={{90,-40},{110,40}})));
    equation
      //the main equation
      rr = kC * ((p * products.u) - (s * substrates.u));

      //reaction molar rates
      rr*s = -substrates.q;
      rr*p = products.q;

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Rectangle(
              extent={{-100,-30},{100,30}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-100,40},{100,66}},
              textString="%name",
              lineColor={0,0,255}),
            Polygon(
              points={{-60,6},{-60,4},{54,4},{54,4},{18,14},{18,6},{-60,6}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{54,-8},{54,-6},{-60,-6},{-60,-6},{-24,-16},{-24,-8},{54,-8}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-170,-70},{-22,-44}},
              lineColor={0,0,0},
            textString="%s"),
            Text(
              extent={{12,-68},{160,-42}},
              lineColor={0,0,0},
            textString="%p")}),
        Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &LT;-&GT; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</b></sub> </p>
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
      extends Interfaces.PartialSubstanceInSolution(final stateOfMatter, final substanceData(
        MolarWeight=5.4857990946e-7,
        z=-1,
        DfH_25degC=0,
        DfG_25degC_1bar=0,
        Cp=0,
        density=1e20));

      Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
          Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
              extent={{-110,-10},{-90,10}})));

    equation
      //electric
      pin.v = electricPotential;
      pin.i + z*Modelica.Constants.F*port_a.q = 0;

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
              lineColor={0,0,255})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ElectronTransfer;

    model Diffusion "Solute diffusion"
      extends Icons.Diffusion;
      extends Interfaces.OnePortParallel;
      extends Interfaces.ConditionalKinetics;

    equation
      port_b.q = kC * (port_b.u - port_a.u);

       annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
    end Diffusion;

    model GasSolubility "Henry's law of gas solubility in liquid."

      extends Icons.GasSolubility;

      parameter Boolean useWaterCorrection = true
      "Are free Gibbs energy of aqueous formation shifted by 10 kJ/mol?"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));

    Interfaces.SubstancePort_b gas_port "Gaseous solution"
      annotation (Placement(transformation(extent={{-10,90},{10,110}})));

    Interfaces.SubstancePort_b liquid_port "Dissolved in liquid solution"
      annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
          iconTransformation(extent={{-10,-110},{10,-90}})));

       extends Interfaces.ConditionalKinetics;

    equation
      gas_port.q + liquid_port.q = 0;

      // the main equation
      liquid_port.q = kC *(liquid_port.u - gas_port.u - (if useWaterCorrection then Modelica.Constants.R*(298.15)*log(0.01801528) else 0));

       annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
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
      extends Interfaces.OnePortParallel;
      extends Interfaces.ConditionalKinetics;

    equation
      //the main equation
      port_a.q = kC * (port_a.u - port_b.u);

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
              lineColor={0,0,255},
            origin={69,2},
            rotation=90)}));
    end Membrane;

    model SubstancePump "Prescribed sunstance molar flow"
      extends Interfaces.OnePortParallel;
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
              lineColor={0,0,255},
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

      parameter Integer NumberOfSubunits=1
      "Number of independent subunits occurring in macromolecule";

    Interfaces.SubstancePort_b subunits[NumberOfSubunits]
      "Subunits of macromolecule" annotation (Placement(transformation(extent={
              {-30,90},{-10,110}}), iconTransformation(extent={{-30,90},{-10,
              110}})));

      Interfaces.SolutionPort solution annotation (Placement(transformation(extent={{-70,
                -110},{-50,-90}}),
            iconTransformation(extent={{-70,-110},{-50,-90}})));

        Modelica.SIunits.AmountOfSubstance nm
      "Amount of the macromolecule (all form in the conformation)";
        Modelica.SIunits.MoleFraction xm
      "Mole fraction of the macromolecule (all form of in the conformation)";

  public
      Interfaces.SolutionPort subunitSolution
      "The port to connect all subunits"
        annotation (Placement(transformation(extent={{-70,92},{-50,112}}),
            iconTransformation(extent={{30,90},{50,110}})));
    Interfaces.SubstancePort_a port_a annotation (Placement(transformation(
            extent={{90,-110},{110,-90}}), iconTransformation(extent={{90,-110},
              {110,-90}})));
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

      //properties from subunits
      subunitSolution.dH + solution.dH = 0;
      subunitSolution.i + solution.i = 0;
      subunitSolution.Qj + solution.Qj = 0;
      subunitSolution.Ij + solution.Ij = 0;

      //properties of macromolecule as a whole
      subunitSolution.nj + solution.nj*NumberOfSubunits = 0;
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
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
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
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics={                                                        Text(
              extent={{-22,-106},{220,-140}},
              lineColor={0,0,255},
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
      SubstancePump substancePump(useSubstanceFlowInput=true)
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
    equation
      product.u1=q;
      product1.u2=q;
      greaterThreshold.u=q;

      connect(port_b, moleFractionSensor1.port_a) annotation (Line(
          points={{-100,0},{-76,0}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(moleFractionSensor.port_a, port_a) annotation (Line(
          points={{76,0},{100,0}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(moleFractionSensor1.solution, solution) annotation (Line(
          points={{-60,-10},{-60,-100}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(solution, moleFractionSensor.solution) annotation (Line(
          points={{-60,-100},{60,-100},{60,-10}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(substancePump.substanceFlow, switch1.y) annotation (Line(
          points={{0,-60},{0,-49},{-2.22045e-015,-49}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(switch1.u2, greaterThreshold.y) annotation (Line(
          points={{2.22045e-015,-26},{0,-26},{0,-15}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(product1.u1, moleFractionSensor.moleFraction) annotation (Line(
          points={{42,-32},{50,-32},{50,0},{56,0}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(product.u2, moleFractionSensor1.moleFraction) annotation (Line(
          points={{-42,-32},{-50,-32},{-50,0},{-56,0}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(port_b, substancePump.port_a) annotation (Line(
          points={{-100,0},{-86,0},{-86,-64},{-14,-64}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(substancePump.port_b, port_a) annotation (Line(
          points={{6,-64},{84,-64},{84,0},{100,0}},
          color={158,66,200},
          smooth=Smooth.None));
      connect(product.y, switch1.u1) annotation (Line(
          points={{-19,-26},{-8,-26}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(product1.y, switch1.u3) annotation (Line(
          points={{19,-26},{8,-26}},
          color={0,0,127},
          smooth=Smooth.None));
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
              lineColor={0,0,255},
              origin={2,-74},
              rotation=180)}),
        Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><h4><font color=\"#008000\">Bidirectional mass flow by concentration</font></h4></p>
<p>Possible field values: </p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0.1\"><tr>
<td></td>
<td><p align=\"center\"><h4>forward flow</h4></p></td>
<td><p align=\"center\"><h4>backward flow</h4></p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>solutionFlow</h4></p></td>
<td><p align=\"center\">&GT;=0</p></td>
<td><p align=\"center\">&LT;=0</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>q_in.q</h4></p></td>
<td><p align=\"center\">=solutionFlow*q_in.conc</p></td>
<td><p align=\"center\">=-q_out.q</p></td>
</tr>
<tr>
<td><p align=\"center\"><h4>q_out.q</h4></p></td>
<td><p align=\"center\">=-q_in.q</p></td>
<td><p align=\"center\">=solutionFlow*q_out.conc</p></td>
</tr>
</table>
<br/>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics));
    end Stream;

    model FluidAdapter

      constant String substanceNames[:] annotation (Evaluate=true);

      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium (substanceNames=substanceNames)
      "Medium model"   annotation (choicesAllMatching=true);                                               //reducedX=true, fixedX=false,

      Modelica.Fluid.Interfaces.FluidPort_a fluid(redeclare package Medium = Medium)
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      Interfaces.SubstancePorts_b substances[n] annotation (Placement(transformation(
              extent={{-110,-40},{-90,40}}), iconTransformation(extent={{-110,-40},{
                -90,40}})));
      Interfaces.SolutionPort solution
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}}),
            iconTransformation(extent={{-10,-60},{10,-40}})));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics), Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={Line(
              points={{-90,0},{90,0}},
              color={158,66,200},
              thickness=1,
              smooth=Smooth.None)}));
      replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby
      Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
         annotation (choicesAllMatching = true);

      parameter stateOfMatter.SubstanceData substanceData[n]
      "Definition of the substance"
         annotation (choicesAllMatching = true);

      Modelica.SIunits.MoleFraction x[n] "Mole fraction of the substance";

      Modelica.SIunits.ActivityOfSolute a[n]
      "Activity of the substance (mole-fraction based)";

      Modelica.SIunits.ActivityCoefficient gamma[n]
      "Activity coefficient of the substance";

      Modelica.SIunits.ChargeNumberOfIon z[n] "Charge number of ion";

  protected
      parameter Integer n=Medium.nS "Number of substances";
      Modelica.SIunits.MolarMass molarMass[n] "Molar mass of the substance";

      Modelica.SIunits.Temperature temperature(start=298.15)
      "Temperature of the solution";

      Modelica.SIunits.Pressure pressure(start=100000)
      "Pressure of the solution";

      Modelica.SIunits.ElectricPotential electricPotential(start=0)
      "Electric potential of the solution";

      Modelica.SIunits.MoleFraction moleFractionBasedIonicStrength(start=0)
      "Ionic strength of the solution";

      Modelica.SIunits.AmountOfSubstance amountOfSolution
      "Amount of all solution particles";

    //  Medium.SpecificHeatCapacity Cp;

    //  Medium.ThermodynamicState thermodynamicState;
    equation
      //fluid
      fluid.m_flow + molarMass*substances.q = 0;
      fluid.p = pressure;
      solution.dH = if
                      (fluid.m_flow >0) then inStream(fluid.h_outflow) * fluid.m_flow else 0;
      fluid.h_outflow = 0; //outgoing substance does change solution.dH using substance components

      for i in 1:n-1 loop
       // actualStream(fluid.Xi_outflow[i]) = molarMass[i] * substances[i].q ./ fluid.m_flow;
       // fluid.Xi_outflow[i] = x[i] * molarMass[i] / (x*molarMass);

        substances[i].q = if
                            (fluid.m_flow > 0) then fluid.m_flow*inStream(fluid.Xi_outflow[i]) / molarMass[i] else fluid.m_flow * (x[i] / (x*molarMass));
      end for;

     // thermodynamicState = Medium.setState_pTX(pressure, temperature, actualStream(fluid.Xi_outflow));
     // Cp = Medium.specificHeatCapacityCp(thermodynamicState);

      //substances
      for i in 1:n loop
       gamma[i] = stateOfMatter.activityCoefficient(substanceData[i],temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
       z[i] = stateOfMatter.chargeNumberOfIon(substanceData[i],temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
       molarMass[i] = stateOfMatter.molarMass(substanceData[i],temperature,pressure,electricPotential,moleFractionBasedIonicStrength);

       a[i] = gamma[i]*x[i];
       substances[i].u = stateOfMatter.chemicalPotentialPure(
        substanceData[i],
        temperature,
        pressure,
        electricPotential,
        moleFractionBasedIonicStrength)
        + Modelica.Constants.R*temperature*log(a[i])
        + z[i]*Modelica.Constants.F*electricPotential;
      end for;

      //solution
      temperature = solution.T;
      pressure = solution.p;
      electricPotential = solution.v;
      amountOfSolution = solution.n;
      moleFractionBasedIonicStrength = solution.I;

      solution.i = 0;
      solution.dV = 0;
      solution.Gj = 0;
      solution.nj = 0;
      solution.mj = 0;
      solution.Qj = 0;
      solution.Ij = 0;
      solution.Vj = 0;
    end FluidAdapter;
  end Components;


  package Sensors "Chemical sensors"
    extends Modelica.Icons.SensorsPackage;

    model MolarFlowSensor "Measure of molar flow"

      extends Modelica.Icons.RotationalSensor;
      extends Interfaces.OnePortSerial;

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
              color={127,0,127},
              smooth=Smooth.None),
            Line(
              points={{70,10},{90,10}},
              color={127,0,127},
              smooth=Smooth.None),
            Line(
              points={{-90,10},{-70,10}},
              color={127,0,127},
              smooth=Smooth.None),
            Line(
              points={{-90,-10},{-70,-10}},
              color={127,0,127},
              smooth=Smooth.None),
            Text(
              extent={{-31,-5},{28,-64}},
              lineColor={0,0,0},
              textString="dn")}));
    end MolarFlowSensor;

    model MoleFractionSensor "Measure of mole fraction"
      extends Modelica.Icons.RotationalSensor;
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
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MoleFractionSensor;

    model ElectroChemicalPotentialSensor
    "Measure of electro-chemical potential"
      extends Modelica.Icons.RotationalSensor;

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
      port_a.q = 0;

      port_a.u = u;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="u"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ElectroChemicalPotentialSensor;

    model MolalitySensor "Measure of molality of the substance"
      extends Modelica.Icons.RotationalSensor;
      extends Interfaces.PartialSubstanceSensor;

      parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionPer1kgOfSolvent = 55.508
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
      constant Modelica.SIunits.Mass KG=1;
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
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MolalitySensor;

    model MolarConcentrationSensor "Measure of molarity of the substance"
      extends Modelica.Icons.RotationalSensor;
      extends Interfaces.PartialSubstanceSensor;

    parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionInOneLiter = 55.508
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
      constant Modelica.SIunits.Volume L=0.001;
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
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MolarConcentrationSensor;

    model MassFractionSensor "Measure of mass fraction of the substance"
      extends Modelica.Icons.RotationalSensor;
      extends Interfaces.PartialSubstanceSensor;

    parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionInOneKilogram = 55.508
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

      x=(massFraction/molarMass) / AmountOfSolutionInOneKilogram;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="mx"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MassFractionSensor;

    model PartialPressureSensor
    "Measure of partial pressure of the substance in gaseous solution"
      extends Modelica.Icons.RotationalSensor;
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
              color={127,0,127},
              smooth=Smooth.None)}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PartialPressureSensor;

    model DissociationCoefficient
    "Meassure dissociation coefficient (mole fraction based) for pure substances"
      extends Modelica.Icons.TranslationalSensor;

      parameter Modelica.SIunits.Temperature T=298.15 "Temperature";

      parameter Modelica.SIunits.AmountOfSubstance n=55.508*m
      "Amount of all substances in solution per one liter of solution";

      parameter Modelica.SIunits.Mass m=0.997
      "Mass of solvent per one liter of solution";

      parameter Integer nS=1 "Number of substrates types"
        annotation ( HideResult=true);

      parameter Modelica.SIunits.StoichiometricNumber s[nS]=ones(nS)
      "Stoichiometric reaction coefficient for substrates"
        annotation (HideResult=true);

      parameter Integer nP=1 "Number of products types"
        annotation ( HideResult=true);

      parameter Modelica.SIunits.StoichiometricNumber p[nP]=ones(nP)
      "Stoichiometric reaction coefficients for products"
        annotation (HideResult=true);

    Interfaces.SubstancePort_b products[nP] "Products"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    Interfaces.SubstancePort_b substrates[nS] "Substrates"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

      Modelica.SIunits.MolarEnergy DrG "Free Gibbs energy of reaction";

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
    equation
      substrates.q = zeros(nS);
      products.q = zeros(nP);

      DrG = ((p * products.u) - (s * substrates.u));

      DissociationCoefficient_MoleFractionBased = exp(-DrG/(Modelica.Constants.R*T));

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
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &LT;-&GT; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</b></sub> </p>
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

  end Sensors;


  package Sources "Chemical sources"
    extends Modelica.Icons.SourcesPackage;

    model PureSubstance "Constant source of pure substance"
      extends Interfaces.PartialSubstance;

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
      "Electric potential";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
      "Ionic strength";

    equation
      x = 1;

      //the solution
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
              lineColor={0,0,255}),
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
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Pressure PartialPressure=0
      "Fixed partial pressure if usePartialPressureInput=false"
        annotation (HideResult=true, Dialog(enable=not usePartialPressureInput));

      parameter Modelica.SIunits.Pressure TotalPressure=101325
      "Total pressure of the whole gaseous solution";

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
      "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
      "Electric potential";

      Modelica.Blocks.Interfaces.RealInput partialPressure(start=
            PartialPressure, final unit="Pa")=p if usePartialPressureInput
      "Partial pressure of gas = total pressure * gas fraction"
        annotation (HideResult=true,Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.Pressure p "Current partial pressure";

      parameter Modelica.SIunits.Volume Volume = 0.001
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
              smooth=Smooth.None,
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
              lineColor={0,0,255}),
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

       parameter Real Molality(final unit="mol/kg") = 1e-8
      "Fixed molality of the substance if useMolalityInput=false"
        annotation (HideResult=true, Dialog(enable=not useMolalityInput));

       parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionPer1KgSolvent = 55.508
      "Amount of all particles in the solution per one kilogram of solvent";

        parameter Boolean useMolalityInput = false
      "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
      "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
      "Electric potential";

      Modelica.Blocks.Interfaces.RealInput molalityInput(start=Molality,final unit="mol/kg")=n/KG if
           useMolalityInput
        annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.AmountOfSubstance n "Current amount of the substance";

  protected
      constant Modelica.SIunits.Mass KG=1;
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
              lineColor={0,0,255}),
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

       parameter Real MolarConcentration(final unit="mol/m3", displayUnit="mol/l") = 1e-8
      "Fixed molarity of the substance if useMolarityInput=false"
        annotation (HideResult=true, Dialog(enable=not useMolarityInput));

       parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 55.508
      "Amount of all particles in the solution one liter of solvent";

        parameter Boolean useMolarityInput = false
      "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

       parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
      "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
      "Electric potential";

      Modelica.Blocks.Interfaces.RealInput molarConcentrationInput(start=MolarConcentration,final unit="mol/m3", displayUnit="mol/l")=n/L if
           useMolarityInput
        annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

      Modelica.SIunits.AmountOfSubstance n "Current amount of the substance";

  protected
      constant Modelica.SIunits.Volume L=0.001;
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
              lineColor={0,0,255}),
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

       parameter Modelica.SIunits.MoleFraction MoleFraction = 1e-8
      "Fixed mole fraction of the substance if useMoleFractionInput=false"
        annotation (HideResult=true, Dialog(enable=not useMoleFractionInput));

        parameter Boolean useMoleFractionInput = false
      "Is mole fraction of the substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Temperature Temperature=298.15 "Temperature";
      parameter Modelica.SIunits.Pressure Pressure=101325 "Pressure";
      parameter Modelica.SIunits.MoleFraction MoleFractionBasedIonicStrength=0
      "Ionic strength";
      parameter Modelica.SIunits.ElectricPotential ElectricPotential=0
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
              lineColor={0,0,255}),
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

       parameter Modelica.SIunits.ChemicalPotential U = 1e-8
      "Fixed electro-chemical potential of the substance if usePotentialInput=false"
        annotation (HideResult=true, Dialog(enable=not usePotentialInput));

       parameter Boolean usePotentialInput = false
      "Is electro-chemical potential of the substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      Modelica.Blocks.Interfaces.RealInput uInput(final unit="J/mol")=port_a.u if
           usePotentialInput annotation (HideResult=true, Placement(transformation(
              extent={{-120,-20},{-80,20}})));

    Interfaces.SubstancePort_a port_a
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));
    equation
       if not usePotentialInput then
         port_a.u=U;
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
              lineColor={0,0,255}),
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

    equation
      port_b.q = - q;

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
              extent={{-82,-82},{90,-58}},
              textString="%name",
              lineColor={0,0,255})}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceInflow;

    model SubstanceOutflow "Molar pump of substance out of system"
      extends Interfaces.ConditionalSubstanceFlow;

    Interfaces.SubstancePort_b port_a "Inflow"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

    equation
      port_a.q = q;

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
              extent={{-82,-82},{90,-58}},
              textString="%name",
              lineColor={0,0,255})}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstanceOutflow;

    model Clearance "Physiological Clearance"
     extends Interfaces.ConditionalSolutionFlow(final SolutionFlow=Clearance/K);
     extends Interfaces.PartialSubstanceSensor;

      parameter Modelica.SIunits.VolumeFlowRate Clearance=0
      "Physiological clearance of the substance if useSolutionFlowInput=false"
        annotation (HideResult=true, Dialog(enable=not useSolutionFlowInput));

      parameter Real K(unit="1")=1
      "Coefficient such that Clearance = K*solutionFlow";

      Modelica.SIunits.MolarFlowRate molarClearance "Current molar clearance";

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
              lineColor={0,0,255}),
            Text(
              extent={{-100,-30},{100,-50}},
              lineColor={0,0,0},
              textString="K=%K")}),        Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Clearance;

    model Degradation "Degradation of the substance"
      extends Interfaces.PartialSubstanceSensor;

      parameter Modelica.SIunits.Time HalfTime
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
              lineColor={0,0,255}),
            Text(
              extent={{-100,54},{100,28}},
              lineColor={0,0,0},
              textString="t1/2 = %HalfTime s"),
            Polygon(
              points={{54,24},{54,-24},{44,-22},{44,22},{54,24}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{30,20},{30,-20},{20,-18},{20,18},{30,20}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{8,16},{8,-16},{-2,-14},{-2,14},{8,16}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-12,12},{-12,-12},{-22,-10},{-22,10},{-12,12}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-34,8},{-34,-8},{-44,-6},{-44,6},{-34,8}},
              lineColor={0,0,127},
              smooth=Smooth.None,
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-56,4},{-56,-4},{-66,-2},{-66,2},{-56,4}},
              lineColor={0,0,127},
              smooth=Smooth.None,
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
<td>2013-2015</td>
</tr>
</table>
</html>"));
    end Degradation;

    model Buffer
    "Source of substance bounded to constant amount of buffer to reach linear dependence between concentration and electrochemical potential"
      extends Icons.Buffer;
           extends Interfaces.PartialSubstanceInSolution(a(start = a_start));

         parameter Modelica.SIunits.MoleFraction a_start=1e-7
      "Initial value of mole fraction of the buffered substance";

         parameter Modelica.SIunits.AmountOfSubstance BufferValue = 0.001
      "Fixed buffer value (slope between amount of buffered substance and -log10(activity)) if useBufferValueInput=false"
          annotation (HideResult=true, Dialog(enable=not useMoleFractionInput));

         parameter Boolean useBufferValueInput = false
      "Is buffer value of the substance an input?"
          annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

          extends Interfaces.ConditionalKinetics;

          Real bufferValue(final unit="1");

        Modelica.Blocks.Interfaces.RealInput bufferValueInput(
          final unit="mol/mol",
          start=BufferValue)=bufferValue if
             useBufferValueInput annotation (HideResult=true, Placement(transformation(
                extent={{-120,-20},{-80,20}})));

          Real xref;
        Modelica.SIunits.AmountOfSubstance nFreeBuffer(start=-log10(a_start)*BufferValue);
        Modelica.SIunits.MoleFraction xFreeBuffer;

  protected
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
      solution.dH = molarEnthalpy*port_a.q - der(molarEnthalpy)*nFreeBuffer;
      solution.i = Modelica.Constants.F * z * port_a.q - Modelica.Constants.F*der(z)*nFreeBuffer;
      solution.dV = molarVolume * port_a.q - der(molarVolume)*nFreeBuffer;

      //extensive properties
      solution.nj=0;
      solution.mj=-nFreeBuffer*molarMass;
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
                lineColor={0,0,255})}),
          Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Buffer;
  end Sources;


  package Interfaces "Chemical interfaces"
    extends Modelica.Icons.InterfacesPackage;

    connector SubstancePort
    "Electro-chemical potential and molar change of the substance in the solution"

      Modelica.SIunits.ChemicalPotential u
      "Electro-chemical potential of the substance in the solution";

      flow Modelica.SIunits.MolarFlowRate q "Molar change of the substance";

      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Definition of electro-chemical potential of the substance:</p>
<p><b>u<sub>e-ch</sub>(x,T,V) = u&deg;(T) + R*T*ln(gamma*x) + z*F*V</b></p>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>V .. eletric potential of the substance</p>
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
              lineThickness=1),
       Text(extent = {{-160,110},{40,50}}, lineColor={172,72,218},   textString = "%name")}),
        Documentation(info="<html>
<p>Chemical port with external definition of the substance outside the component.</p>
</html>",
        revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end SubstancePort_b;

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

    SubstancePort_a port_a "The substance"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      replaceable package stateOfMatter = Incompressible                    constrainedby
      StateOfMatter
      "Substance model to translate data into substance properties"
         annotation (choicesAllMatching = true);

      parameter stateOfMatter.SubstanceData substanceData
      "Definition of the substance"
         annotation (choicesAllMatching = true);

      Modelica.SIunits.MoleFraction x "Mole fraction of the substance";

      Modelica.SIunits.ActivityOfSolute a
      "Activity of the substance (mole-fraction based)";

  protected
      Modelica.SIunits.ActivityCoefficient gamma
      "Activity coefficient of the substance";

      Modelica.SIunits.ChargeNumberOfIon z "Charge number of ion";

      Modelica.SIunits.Temperature temperature(start=298.15)
      "Temperature of the solution";

      Modelica.SIunits.Pressure pressure(start=100000)
      "Pressure of the solution";

      Modelica.SIunits.ElectricPotential electricPotential(start=0)
      "Electric potential of the solution";

      Modelica.SIunits.MoleFraction moleFractionBasedIonicStrength(start=0)
      "Ionic strength of the solution";

      Modelica.SIunits.MolarMass molarMass "Molar mass of the substance";

      Modelica.SIunits.MolarEnthalpy molarEnthalpy
      "Molar enthalpy of the substance";

      Modelica.SIunits.MolarEntropy molarEntropyPure
      "Molar entropy of the pure substance";

      Modelica.SIunits.ChemicalPotential u0
      "Chemical potential of the pure substance";

      Modelica.SIunits.ChemicalPotential uPure
      "Electro-Chemical potential of the pure substance";

      Modelica.SIunits.MolarVolume molarVolume "Molar volume of the substance";

      Modelica.SIunits.MolarVolume molarVolumePure
      "Molar volume of the pure substance";

      Modelica.SIunits.MolarVolume molarVolumeExcess
      "Molar volume excess of the substance in solution (typically it is negative as can be negative)";

    //  Modelica.SIunits.MolarHeatCapacity molarHeatCapacityCp
    //    "Molar heat capacity of the substance at constant pressure";

    equation
      //aliases
      gamma = stateOfMatter.activityCoefficient(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
      z = stateOfMatter.chargeNumberOfIon(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
      molarMass = stateOfMatter.molarMass(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);

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

      annotation (
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        experiment);
    end PartialSubstance;

    partial model PartialSubstanceInSolution
    "Substance properties for components, where the substance is connected with the solution"

      SolutionPort            solution
      "To connect substance with solution, where is pressented"                                  annotation (Placement(transformation(
              extent={{-70,-110},{-50,-90}}),iconTransformation(extent={{-70,-110},{
                -50,-90}})));

      extends PartialSubstance;

  protected
       Modelica.SIunits.AmountOfSubstance amountOfSolution
      "Amount of all solution particles";

    equation
      //aliases
      temperature = solution.T;
      pressure = solution.p;
      electricPotential = solution.v;
      amountOfSolution = solution.n;
      moleFractionBasedIonicStrength = solution.I;

    end PartialSubstanceInSolution;

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

     replaceable record SubstanceData
      "Definition data of the chemical substance"
     end SubstanceData;

     replaceable function activityCoefficient
      "Return activity coefficient of the substance in the solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Real activityCoefficient "Activity Coefficient";
     end activityCoefficient;

     replaceable function chargeNumberOfIon
      "Return charge number of the substance in the solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.ChargeNumberOfIon chargeNumberOfIon
        "Charge number of ion";
     end chargeNumberOfIon;

     replaceable function molarEnthalpyElectroneutral
      "Molar enthalpy of the substance in electroneutral solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarEnthalpy molarEnthalpyElectroneutral
        "Molar enthalpy";
     end molarEnthalpyElectroneutral;

     function molarEnthalpy
      "Molar enthalpy of the substance with electric potential dependence"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarEnthalpy molarEnthalpy "Molar enthalpy";
     algorithm
        molarEnthalpy := molarEnthalpyElectroneutral(substanceData,T,p) +
             Modelica.Constants.F*chargeNumberOfIon(substanceData,T,p,v,I)*v;
     end molarEnthalpy;

     replaceable function molarEntropyPure
      "Molar entropy of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarEntropy molarEntropyPure
        "Molar entropy of the pure substance";
     end molarEntropyPure;

      function molarEntropy "Molar entropy of the substance in the solution"
            extends Modelica.Icons.Function;
            input Modelica.SIunits.ChemicalPotential u
        "Electro-chemical potential of the substance";
            input SubstanceData substanceData "Data record of substance";
            input Modelica.SIunits.Temperature T=298.15 "Temperature";
            input Modelica.SIunits.Pressure p=100000 "Pressure";
            input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
            input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
            output Modelica.SIunits.MolarEntropy molarEntropy "Molar entropy";
      algorithm
          molarEntropy :=  (u - molarEnthalpy(substanceData,T,p,v,I))/T;
      end molarEntropy;

     function chemicalPotentialPure "Chemical potential of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.ChemicalPotential chemicalPotentialPure
        "Base chemical potential";
     algorithm
         chemicalPotentialPure :=  molarEnthalpyElectroneutral(substanceData,T,p) - T*molarEntropyPure(substanceData,T,p,v,I);
     end chemicalPotentialPure;

     function electroChemicalPotentialPure
      "Electro-chemical potential of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.ChemicalPotential electroChemicalPotentialPure
        "Base electro-chemical potential";
     algorithm
      electroChemicalPotentialPure := chemicalPotentialPure(
           substanceData,
           T,
           p,
           v,
           I) + Modelica.Constants.F*chargeNumberOfIon(substanceData,T,p,v,I)*v;
     end electroChemicalPotentialPure;

     replaceable function molarMass "Molar mass of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarMass molarMass "Molar mass";
     end molarMass;

     replaceable function molarVolumePure "Molar volume of the pure substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarVolume molarVolumePure "Molar volume";
     end molarVolumePure;

     function molarVolumeExcess
      "Excess molar volume of the substance in the solution"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarVolume molarVolumeExcess
        "Excess molar volume of the substance in the solution";
     algorithm
        molarVolumeExcess := molarVolumePure(substanceData,T,p,v,I)*
           log(activityCoefficient(substanceData,T,p,v,I)); //zero if activityCoefficient==1
     end molarVolumeExcess;

     replaceable function molarVolume "Molar volume of the substance"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarVolume molarVolume "Molar volume";
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
     end molarVolume;

      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end StateOfMatter;

    package IdealGas "Ideal gas as basic state of matter"
       extends StateOfMatter;

       redeclare replaceable record extends SubstanceData "Base substance data"

          parameter Modelica.SIunits.MolarMass MolarWeight(displayUnit="kDa")=0.01801528
        "Molar weight of the substance in kg/mol or kDa";

          parameter Modelica.SIunits.ChargeNumberOfIon z=0
        "Charge number of the substance (e.g. 0..uncharged, -1..electron, +2..Ca^2+)";

          parameter Modelica.SIunits.MolarEnergy DfH_25degC(displayUnit="kJ/mol")=0
        "Enthalpy of formation of the substance at 25°C";

          parameter Modelica.SIunits.MolarEnergy DfG_25degC_1bar(displayUnit="kJ/mol")=0
        "Gibbs enerfy of formation of the substance at 25°C,1bar";

          parameter Modelica.SIunits.ActivityCoefficient gamma=1
        "Activity coefficient of the substance";

          parameter Modelica.SIunits.MolarHeatCapacity Cp = 33.6
        "Molar heat capacity of the substance";

          parameter String References[:]={""}
        "References of these thermodynamical values";

        annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
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
         molarEnthalpyElectroneutral := substanceData.DfH_25degC
           +(T-298.15)*(substanceData.Cp+Modelica.Constants.R);
     end molarEnthalpyElectroneutral;

     redeclare function extends molarEntropyPure
      "Molar entropy of the pure substance"
     algorithm
       //molarEntropyPure := ((substanceData.DfH - substanceData.DfG_25degC_1bar)/298.15)
       //+ (substanceData.Cp+Modelica.Constants.R)*log(T/298.15);

         //Molar entropy:
         // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
         // - pressure shift: to reach the ideal gas equation at constant temperature Vm*dP = -T*dS (small amount of work)
         molarEntropyPure := (substanceData.Cp+Modelica.Constants.R)*log(T/298.15) - Modelica.Constants.R*log(p/100000) + ((substanceData.DfH_25degC
          - substanceData.DfG_25degC_1bar)/298.15);

         //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
         //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
     end molarEntropyPure;

     redeclare function extends molarMass "Molar mass of the substance"
     algorithm
         molarMass := substanceData.MolarWeight;
     end molarMass;

     redeclare function extends molarVolumePure
      "Molar volume of the pure substance"
     algorithm
         molarVolumePure := Modelica.Constants.R*T/p; //ideal gas
     end molarVolumePure;

    /* redeclare function extends molarHeatCapacityCp
    "Molar heat capacity of the substance at constant pressure"
 algorithm
     molarHeatCapacity := substanceData.Cp;
     end molarHeatCapacityCp; */

      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end IdealGas;

    package Incompressible "Incompressible as basic state of matter"
       extends StateOfMatter;

       redeclare replaceable record extends SubstanceData "Base substance data"

          parameter Modelica.SIunits.MolarMass MolarWeight(displayUnit="kDa")=0.01801528
        "Molar weight of the substance in kg/mol or kDa";

          parameter Modelica.SIunits.ChargeNumberOfIon z=0
        "Charge number of the substance (e.g. 0..uncharged, -1..electron, +2..Ca^2+)";

          parameter Modelica.SIunits.MolarEnergy DfH_25degC(displayUnit="kJ/mol")=0
        "Enthalpy of formation of the substance at 25°C";

          parameter Modelica.SIunits.MolarEnergy DfG_25degC_1bar(displayUnit="kJ/mol")=0
        "Gibbs enerfy of formation of the substance at 25°C,1bar";

          parameter Modelica.SIunits.ActivityCoefficient gamma=1
        "Activity coefficient of the substance";

          parameter Modelica.SIunits.MolarHeatCapacity Cp = 75.3
        "Molar heat capacity of the substance at constant pressure";

        //      parameter Modelica.SIunits.MolarHeatCapacity Cv = Cp
        //      "Molar heat capacity of the substance at constant volume";

          parameter Modelica.SIunits.Density density(displayUnit="kg/dm3")=997.0479
        "Density of the pure substance (default density of water at 25degC)";

          parameter String References[:]={""}
        "References of these thermodynamical values";

        annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
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

         molarEnthalpyElectroneutral :=  substanceData.DfH_25degC
         + (T - 298.15) * substanceData.Cp;
      //   - (p - 100000) * molarVolumePure(substanceData,T,p,v,I);
     end molarEnthalpyElectroneutral;

      redeclare function extends molarEntropyPure
      "Molar entropy of the pure substance"
      algorithm
         //molarEntropyPure := ((substanceData.DfH - substanceData.DfG_25degC_1bar)/298.15)
         //+ substanceData.Cv*log(T/298.15);

         //Molar entropy shift:
         // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
         // - pressure shift: with constant molar volume at constant temperature Vm*dP = -T*dS (small amount of work)
         molarEntropyPure := substanceData.Cp*log(T/298.15) - (molarVolumePure(
           substanceData,
           T,
           p,
           v,
           I)/T)*(p - 100000) + ((substanceData.DfH_25degC - substanceData.DfG_25degC_1bar)/298.15);

         //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
         //As data: http://www1.lsbu.ac.uk/water/water_phase_diagram.html
         //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
      end molarEntropyPure;

     redeclare function extends molarMass "Molar mass of the substance"
     algorithm
         molarMass := substanceData.MolarWeight;
     end molarMass;

     redeclare function extends molarVolumePure
      "Molar volume of the pure substance"
     algorithm
         molarVolumePure := substanceData.MolarWeight/substanceData.density; //incompressible
     end molarVolumePure;

    /* redeclare function extends molarHeatCapacityCp
    "Molar heat capacity of the substance at constant pressure"
 algorithm
     molarHeatCapacity := substanceData.Cp;
 end molarHeatCapacityCp;

 function molarHeatCapacityCv
    "Molar heat capacity of the substance at constant volume"
    extends Modelica.Icons.Function;
   input SubstanceData substanceData "Data record of substance";
   input Modelica.SIunits.Temperature T=298.15 "Temperature";
   input Modelica.SIunits.Pressure p=101325 "Pressure";
   input Modelica.SIunits.ElectricPotential v=0
      "Electric potential of the substance";
   input Modelica.SIunits.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
   output Modelica.SIunits.MolarHeatCapacity molarHeatCapacity
      "Molar heat capacity";
 algorithm
     molarHeatCapacity := substanceData.Cv;
 end molarHeatCapacityCv;
*/
      annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Incompressible;

    connector SolutionPort
    "Only for connecting the one solution their substances. Please, do not use it in different way."

      //enthalpy
      Modelica.SIunits.Temperature T "Temperature of the solution";
      flow Modelica.SIunits.EnthalpyFlowRate dH
      "Internal enthalpy change of the solution";

      //pressure
      Modelica.SIunits.Pressure p "Pressure of the solution";
      flow Modelica.SIunits.VolumeFlowRate dV "Volume change of the solution";

      //electric port
      Modelica.SIunits.ElectricPotential v "Electric potential in the solution";
      flow Modelica.SIunits.ElectricCurrent i "Change of electric charge";

      //Extensive properties of the solution:

      // The extensive quantities here have not the real physical flows.
      // They hack the Kirchhof's flow equation to be counted as the sum from all connected substances in the solution.

      //amount of substances
      Modelica.SIunits.AmountOfSubstance n "Amount of the solution";
      flow Modelica.SIunits.AmountOfSubstance nj
      "Amount of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

      //mass of substances
      Modelica.SIunits.Mass m "Mass of the solution";
      flow Modelica.SIunits.Mass mj
      "Mass of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

      //volume of substances
      Modelica.SIunits.Volume V "Volume of the solution";
      flow Modelica.SIunits.Volume Vj
      "Volume of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

      //free Gibbs energy of substances
      Modelica.SIunits.Energy G "Free Gibbs energy of the solution";
      flow Modelica.SIunits.Energy Gj
      "Free Gibbs energy of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

      //electric charge of the substance
      Modelica.SIunits.ElectricCharge Q "Electric charge of the solution";
      flow Modelica.SIunits.ElectricCharge Qj
      "Electric charge of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

      //ionic strength of substances
      Modelica.SIunits.MoleFraction I
      "Mole fraction based ionic strength of the solution";
      flow Modelica.SIunits.MoleFraction Ij
      "Mole-fraction based ionic strength of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

      annotation (
      defaultComponentName="solution",
      Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Solution port integrates all substances of the solution:</p>
<p>Such as if there are connected together with electric port, thermal port and with port composed with the amont of substance and molar change of substance.</p>
</html>"), Icon(graphics={            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,128,255},
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid)}),
      Diagram(graphics={
       Text(extent={{-160,110},{40,50}},   lineColor={0,128,255},    textString = "%name",
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid),
                      Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={0,128,255},
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid,
              lineThickness=1)}));
    end SolutionPort;

    partial model PartialSolution
    "Chemical solution as homogenous mixture of the substances (only pressure and electric potential are not defined)"

      Modelica.SIunits.Temperature T(start=298.15) "Temperature";

      Modelica.SIunits.Pressure p(start=100000) "Pressure";

      Modelica.SIunits.Volume volume "Current volume of the solution";

      Interfaces.SolutionPort solution "Solution nonflows and flows"
                                      annotation (Placement(
            transformation(extent={{50,-90},{70,-70}}),  iconTransformation(extent={{58,-100},
              {62,-96}})));

  protected
      Modelica.SIunits.Energy freeInternalEnergy(start=0)
      "Free Internal energy of the solution relative to start of the simulation";

      Modelica.SIunits.Entropy freeEntropy
      "Free entropy of the solution relative to start of the simulation";

      Modelica.SIunits.Energy freeEnthalpy
      "Free enthalpy of the solution relative to start of the simulation";

      Modelica.SIunits.Energy freeGibbsEnergy
      "Free Gibbs energy of the solution relative to start of the simulation";

      Modelica.SIunits.HeatFlowRate heatFromEnvironment
      "External heat flow rate";
       Modelica.SIunits.Power workFromEnvironment "External working power";

       Modelica.SIunits.ElectricCharge charge
      "Current electric charge of the solution";
                                              //(start=electricCharge_start)

    initial equation

      freeInternalEnergy = 0;
    equation
      //internal energy
      der(freeInternalEnergy) = heatFromEnvironment + workFromEnvironment;

      heatFromEnvironment + workFromEnvironment = (-solution.dH) - solution.p*(-solution.dV) - volume*der(solution.p);
      //It is the same as: der(freeEnthalpy)=-solution.dH;

      //thermodinamics equations:
      freeInternalEnergy = freeEnthalpy - volume*p; // H=U+p*V
      freeGibbsEnergy = freeEnthalpy - T*freeEntropy; // G=H-T*S

      //aliases
      solution.p = p;
      solution.T = T;
      solution.G = freeGibbsEnergy;
      solution.Q = charge;
      solution.V = volume;

      //Extensive properties of the solution:

      // The extensive quantities here have not the real physical flows.
      // They hack the Kirchhof's flow equation to be counted as the sum from all connected substances in the solution.

      //amount of substances
      solution.n + solution.nj = 0; //total amount of solution is the sum of amounts of each substance

      //mass of substances
      solution.m + solution.mj = 0; //total mass of solution is the sum masses of each substance

      //free Gibs energy
      solution.G + solution.Gj = 0; //total free Gibbs energy of solution is the sum of free Gibbs energies of each substance

      //ionic strength (mole fraction based)
      solution.I + solution.Ij = 0; //total ionic strength of solution is the ionic strengths of each substance

      //electric charge
      solution.Q + solution.Qj = 0; //total electric charge of solution is the sum of charges of each substance

      //volume
      volume + solution.Vj = 0; //total volume of solution is the sum of volumes of each substance

                                                                                                        annotation (
        Documentation(revisions="<html>
<p>2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<h4>amountOfSubstances = &int; MolarFlows</h4>
<h4>mass = &int; massChanges</h4>
<h4>volume = &int; volumeChanges</h4>
<h4>freeEnthalpy = &int; EnthalpyChanges</h4>
<h4>freeEntropy = &int; EntropyChanges</h4>
<h4>freeGibbsEnergy = &int; GibbsEnergyChanges</h4>
<p>Integration of all substances together into one homogenous mixture - the solution.</p>
</html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));
    end PartialSolution;

    partial model OnePortParallel
    "Partial molar flow between two substance definitions"

    SubstancePort_b port_a annotation (Placement(transformation(extent={{-110,-10},
              {-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));
    SubstancePort_b port_b annotation (Placement(transformation(extent={{90,-10},
              {110,10}}), iconTransformation(extent={{90,-10},{110,10}})));
    equation
      port_a.q + port_b.q = 0;

    end OnePortParallel;

    partial model OnePortSerial
    "Partial transfer of substance from substance definition component to another transfer component (such as MolarFlowSensor)"

    SubstancePort_b port_a annotation (Placement(transformation(extent={{-110,-10},
              {-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));
    SubstancePort_a port_b annotation (Placement(transformation(extent={{90,-10},
              {110,10}}), iconTransformation(extent={{90,-10},{110,10}})));
    equation
      port_a.q + port_b.q = 0;

    end OnePortSerial;

    partial model ConditionalSolutionFlow
    "Input of solution molar flow vs. parametric solution molar flow"

      parameter Boolean useSolutionFlowInput = false
      "Is solution flow an input?"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.VolumeFlowRate SolutionFlow=0
      "Volume flow rate of the solution if useSolutionFlowInput=false"   annotation (
          HideResult=true, Dialog(enable=not useSolutionFlowInput));

      parameter Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L=55.508
      "The amount of all particles in one liter of the solution";

      Modelica.Blocks.Interfaces.RealInput solutionFlow(start=SolutionFlow, final unit="m3/s")=
         q*OneLiter/AmountOfSolutionIn1L if useSolutionFlowInput
         annotation ( HideResult=true, Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,40}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,40})));

      Modelica.SIunits.MolarFlowRate q "Current molar solution flow";

  protected
     constant Modelica.SIunits.Volume OneLiter=0.001 "One liter";

    equation
      if not useSolutionFlowInput then
        q*OneLiter/AmountOfSolutionIn1L = SolutionFlow;
      end if;

    end ConditionalSolutionFlow;

    partial model ConditionalSubstanceFlow
    "Input of substance molar flow vs. parametric substance molar flow"

      parameter Boolean useSubstanceFlowInput = false
      "Is substance flow an input?"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.MolarFlowRate SubstanceFlow=0
      "Volumetric flow of Substance if useSubstanceFlowInput=false"   annotation (
          HideResult=true, Dialog(enable=not useSubstanceFlowInput));

      Modelica.Blocks.Interfaces.RealInput substanceFlow(start=SubstanceFlow, final unit="mol/s")=q if
           useSubstanceFlowInput
           annotation (HideResult=true,
           Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={40,40})));

      Modelica.SIunits.MolarFlowRate q "Current Substance flow";
    equation
      if not useSubstanceFlowInput then
        q = SubstanceFlow;
      end if;

    end ConditionalSubstanceFlow;

    partial model ConditionalKinetics
    "Input of kinetics coefficient vs. parametric kinetics coefficient"

      parameter Boolean useKineticsInput = false
      "Is kinetics coefficient as an input?"
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true),Dialog(group="Conditional inputs"));

      parameter Real KC(final unit="mol2.s-1.J-1")=1
      "Chemical kinetics coefficient if useKineticsInput=false"   annotation (
          HideResult=true, Dialog(enable=not useKineticsInput));

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
  end Interfaces;


  package Icons "Icons for chemical models"
    //extends Modelica.Icons.IconsPackage;
    extends Modelica.Icons.Package;

    partial class Diffusion

      annotation (Icon(graphics={Bitmap(extent={{-100,100},{100,-100}}, fileName=
                  "modelica://Chemical/Resources/Icons/diffusion.png")}));

    end Diffusion;

    class Substance

        annotation ( Icon(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics={Bitmap(extent={{-100,100},{100,-100}}, fileName=
                  "modelica://Chemical/Resources/Icons/Substance.png")}));
    end Substance;

    class Speciation

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics={Bitmap(extent={{-100,100},{100,-100}}, fileName=
                  "modelica://Chemical/Resources/Icons/Speciation.png")}));
    end Speciation;

    class GasSolubility

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Bitmap(extent={{-100,100},{100,-100}},
                fileName=
                  "modelica://Chemical/Resources/Icons/GasSolubility.png")}));
    end GasSolubility;

    class Membrane

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Bitmap(extent={{-100,100},{100,-100}},
                fileName="modelica://Chemical/Resources/Icons/membrane.png")}));
    end Membrane;

    class EnzymeKinetics

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Bitmap(extent={{-80,84},{86,-26}},
                fileName=
                  "modelica://Chemical/Resources/Icons/EnzymeKinetics.png")}));
    end EnzymeKinetics;

    class Solution

        annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
          Line(
            points={{-98,90},{-94,96},{-84,98},{84,98},{96,96},{100,92},{98,86},
                {94,80},{94,80},{94,-92},{94,-92},{94,-96},{92,-100},{88,-100},
                {84,-100},{-84,-100},{-88,-100},{-92,-100},{-94,-96},{-94,-92},
                {-94,24},{-94,78},{-94,80},{-98,90}},
            color={127,0,127},
            smooth=Smooth.Bezier,
            pattern=LinePattern.Dot,
            thickness=0.5)}));
    end Solution;

    class Buffer

        annotation ( Icon(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics={Bitmap(extent={{-100,100},{100,-100}}, fileName=
                  "modelica://Chemical/Resources/Icons/buffer.png")}));
    end Buffer;

    class ElectronTransfer

        annotation ( Icon(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics={Bitmap(extent={{-100,100},{100,-100}}, fileName=
                "modelica://Chemical/Resources/Icons/electron.png")}));
    end ElectronTransfer;
    annotation (Documentation(revisions=""));
  end Icons;


  annotation (
preferredView="info",
version="1.1.0",
versionBuild=1,
versionDate="2015-05-20",
dateModified = "2015-05-20 17:14:41Z",
conversion(
  from(version="1.1.0alpha", script="modelica://Chemical/Resources/Scripts/Dymola/Chemical_from_1.0_to_1.1.mos"),
  from(version="1.0.0", script="modelica://Chemical/Resources/Scripts/Dymola/Chemical_from_1.0_to_1.1.mos")),
uses(Modelica(version="3.2.1")),
  Documentation(revisions="<html>
<p>Licensed by Marek Matejak under the Modelica License 2</p>
<p>Copyright &copy; 2008-2015, Marek Matejak, Charles University in Prague.</p>
<p><br><i>This Modelica package is&nbsp;<u>free</u>&nbsp;software and the use is completely at&nbsp;<u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see&nbsp;<a href=\"modelica://Chemical.UsersGuide.ModelicaLicense2\">UsersGuide.ModelicaLicense2</a>&nbsp;or visit&nbsp;<a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>", info="<html>
<p>At firs there was an equilibration of concentrations, but it does not work at all. In reality for almost all electro-chemical processes is equilibrated always the <a href=\"modelica://Chemical.Components.Substance\">electro-chemical potential</a>, not only the concentration.</p>
<p>The pattern is so strong, that the equilibriation of electro-chemical potential can be aplicated for almost all components: chemical reactions, gas dissolution, diffusion, membrane transports, osmotic fluxes, electrochemical cells, electrodes, ..</p>
<p>Please see the <a href=\"modelica://Chemical.UsersGuide.Overview\">overview</a>.</p>
</html>"));

end Chemical;
