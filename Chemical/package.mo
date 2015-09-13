within ;
package Chemical "Library of Electro-Chemical models (version 1.1.0)"
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
<td valign=\"middle\"><p><a href=\"modelica://Chemical.Components.Substance\">Chemical substance</a></p><p>The chemical substance integrates the amount of the chemical substance and from the properties of the connected solution it presents the electro-chemical potential of the substance using the <a href=\"modelica://Chemical.Interfaces.SubstancePort\">SubstancePort</a>.</p><p>There are two basic <a href=\"modelica://Chemical.Interfaces.StateOfMatter\">states of matter</a>: <a href=\"modelica://Chemical.Interfaces.IdealGas\">ideal gas</a> and <a href=\"modelica://Chemical.Interfaces.Incompressible\">incompressible</a> substance. However, the user can easily (re)define their own state of matter by inserting the correct expressions for the pure substance <a href=\"modelica://Chemical.Interfaces.StateOfMatter.activityCoefficient\">activity coefficient</a>, <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarVolumePure\">molar volume</a>, <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarEntropyPure\">molar entropy</a> and <a href=\"modelica://Chemical.Interfaces.StateOfMatter.molarEnthalpyElectroneutral\">molar enthalpy</a>, based on the current solution state (temperature, pressure, electric potential and ionic strength) and the <a href=\"modelica://Chemical.Interfaces.StateOfMatter.SubstanceData\">substance data</a>. The object-oriented design allows users to define the substance data record as part of the state of matter package. Users can select substance parameters according to the state of matter, redefining the getter functions of substance properties.</p><p>The examples work with ideal gases in case of all gaseous substance and incompressible state of matter in case of liquid or solid. The definition data are the molar mass of the substance, the number of charges of the substance, the molar heat capacity of the substance at a constant pressure, free formation enthalpy, free formation Gibbs energy and density (if incompressible) &mdash; all at a temperature of 25&deg;C and pressure 1 bar. Since these parameters are usually recorded in chemical tables at this standard conditions. In this manner, more than 35 real chemical <a href=\"modelica://Chemical.Examples.Substances\">substances</a> in the example package of this chemical library have already been defined. The usage of these predefined substances&rsquo; data is very simple. In the parameter dialog of the chemical substance, the correct record with this data can be selected, as shown in Figure 1.</p><p>This setting is typically the most important setting of each chemical model. All equilibrium coefficients, standard voltages, dissolution coefficients, saturated vapor pressures and so on, are automatically solved using these substance data. As a result, for example, the chemical reaction component only needs to define the stoichiometry coefficients, and the connected substances reach equilibrium at the correct equilibrium coefficient.</p></td>
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
<p>All files in this directory (Physiolibrary) and in all subdirectories, especially all files that build package &QUOT;Physiolibrary&QUOT; are licensed by <u><b>Marek Matejak</b></u> under the <a href=\"http://opensource.org/licenses/BSD-3-Clause\">BSD 3-Clause License</a> (with exception of files &QUOT;Resources/*&QUOT;). </p>
<h4>Licensor:</h4>
<p>Marek Matej&aacute;k,</p>
<p>Hviezdoslavova 632/41,</p>
<p>916 01 Star&aacute; Tur&aacute;, </p>
<p>Slovak Republic,</p>
<p>email: marek@matfyz.cz</p>
<h4><span style=\"color:#008000\">Organization: </span></h4>
<p>Institute of Pathological Physiology, First Faculty of Medicine, Charles University in Prague,</p>
<p>U Nemocnice 5, 128 53 Prague 2, Czech Republic</p>
<p><br><h4>Copyright notices of the files:</h4></p>
<p>Copyright (c) 2008-2015, Marek Matej&aacute;k, Charles University in Prague</p>
<p><br>All rights reserved. </p>
<p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: </p>
<p>1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. </p>
<p>2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. </p>
<p>3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. </p>
<p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS &QUOT;AS IS&QUOT; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
</html>"));
    end License;

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
<li>Draft a new <a href=\"
https://github.com/xogeny/impact/blob/master/resources/docs/modelica2015/paper/impact.md#impact-on-library-developers\">release from &QUOT;master&QUOT;</a> branch with number &QUOT;vX.Y.Z&QUOT; and with release notes. </li>
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

      extends Interfaces.PartialSolution(T(start=temperature_start),p(start=BasePressure));

      parameter Boolean useElectricPort = false "Is electric port pressent?"
      annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

      parameter Boolean ElectricGround = true
      "Is the solution electric potential equal to zero during simulation (if not useElectricPort)?"
        annotation (HideResult=true, Dialog(enable=not useElectricPort));

      parameter Boolean useMechanicPorts = false "Are mechanic ports pressent?"
      annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.SIunits.Area SurfaceArea=0.01
      "Area for surfacePort to connect MultiBody components"
        annotation (HideResult=true, Dialog(enable=useMechanicPorts));

      parameter Boolean isPistonPositionAbsolute=false
      "Relavite position has zero at initial state without force"
        annotation (HideResult=true, Dialog(enable=useMechanicPorts));

      parameter Boolean useThermalPort = false "Is thermal port pressent?"
      annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

      parameter Boolean ConstantTemperature = true
      "Has the solution constant temperature during simulation (if not useThermalPort)?"
         annotation (HideResult=true, Dialog(enable=not useThermalPort));

      parameter Modelica.SIunits.Pressure BasePressure=100000
      "Ambient pressure if useMechanicPort, start pressure or absolute pressure if ConstantPressure"
        annotation (HideResult=true);

      parameter Modelica.SIunits.Temperature temperature_start=298.15
      "Initial temperature of the solution"
         annotation (Dialog(group="Initialization"));

      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort(T=T,Q_flow=heatFromEnvironment) if useThermalPort annotation (
          Placement(transformation(extent={{-70,-90},{-50,-70}}),iconTransformation(
              extent={{-62,-104},{-58,-100}})));
      Modelica.Electrical.Analog.Interfaces.PositivePin electricPin(v=solution.v,i=solution.i) if useElectricPort annotation (Placement(
            transformation(extent={{-70,70},{-50,90}}),    iconTransformation(
              extent={{-62,98},{-58,102}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_a surfaceFlange(f=f,s=top_s) if useMechanicPorts
      "The pressure of solution generate force on prescribed surface."
        annotation (Placement(transformation(extent={{-10,70},{10,90}}),
            iconTransformation(extent={{-2,98},{2,102}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_b bottom(f=-f,s=top_s - ds) if useMechanicPorts
      "Fix of the cilinder on bottom."   annotation (Placement(transformation(
              extent={{-10,-90},{10,-70}}), iconTransformation(extent={{-2,-104},{2,
                -100}})));

  protected
      parameter Modelica.SIunits.Position positionShift(fixed=false)
      "=0 absolute, otherwise negative";
       Modelica.SIunits.Position top_s,ds;
       Modelica.SIunits.Force f;

    initial equation
      T=temperature_start;
      positionShift= if
                       (isPistonPositionAbsolute) then 0 else volume/SurfaceArea;
      //s=volume/SurfaceArea - positionShift;
    equation

      //hydraulic
      ds = volume/SurfaceArea - positionShift;
      workFromEnvironment = -der(f*ds); //=der( (p-p0) * volume)
      solution.p = BasePressure - f/SurfaceArea;
      if not useMechanicPorts then
        f=0;
        top_s=ds; //equivalent for bottom_s==0
      end if;

      //electric
      if (not useElectricPort) and
                                  ElectricGround then
        //Solution connected to ground has zero voltage. However, electric current from the solution can varies.
          solution.v = 0;
      end if;
      if (not useElectricPort) and
                                 (not ElectricGround) then
        //Electrically isolated solution has not any electric current from/to the solution. However, electric potential can varies.
          solution.i = 0;
      end if;

      //thermal
      if (not useThermalPort) and ConstantTemperature then
        //Ideal thermal exchange between environment and solution to reach constant temperature
        der(T)=0;
      end if;
      if (not useThermalPort) and (not ConstantTemperature) then
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
<h4>amountOfSolution = &sum; amountOfSubstances</h4>
<h4>mass = &sum; massOfSubstances</h4>
<h4>volume = &sum; volumeOfSubstances</h4>
<h4>freeGibbsEnergy = &sum; freeGibbsEnergiesOfSubstances</h4>
<p>To calculate the sum of extensive substance's properties is misused the Modelica \"flow\" prefix even there are not real physical flows. </p>
</html>"));
    end Solution;
    extends Modelica.Icons.Package;

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
          transformation(extent={{90,-40},{110,40}}), iconTransformation(extent=
             {{90,-40},{110,40}})));

      parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  protected
      Modelica.SIunits.ChemicalPotential du;
    equation
      //the main equation
      du = ((p * products.u) - (s * substrates.u));
      rr = kC * du * exp(-kE*abs(du));

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
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{54,-8},{54,-6},{-60,-6},{-60,-6},{-24,-16},{-24,-8},{54,-8}},
              lineColor={0,0,0},
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

      parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  protected
      Modelica.SIunits.ChemicalPotential du;
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

      parameter Boolean useWaterCorrection = true
      "Are free Gibbs energy of aqueous formation shifted by 10 kJ/mol?"
      annotation(Evaluate=true, HideResult=true, choices(checkbox=true));

    Interfaces.SubstancePort_b gas_port "Gaseous solution"
      annotation (Placement(transformation(extent={{-10,90},{10,110}})));

    Interfaces.SubstancePort_b liquid_port "Dissolved in liquid solution"
      annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
          iconTransformation(extent={{-10,-110},{10,-90}})));

       extends Interfaces.ConditionalKinetics;

      parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  protected
      Modelica.SIunits.ChemicalPotential du;
    equation
      gas_port.q + liquid_port.q = 0;

      du = (liquid_port.u - gas_port.u - (if useWaterCorrection then Modelica.Constants.R*(298.15)*log(0.01801528) else 0));
      // the main equation
      liquid_port.q = kC * du * exp(-kE*abs(du));

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

      parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  protected
      Modelica.SIunits.ChemicalPotential du;
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
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
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
</html>"));
    end Stream;

    model FluidAdapter
    "Adapter between chemical substances of one homogenous chemical solution and Modelica.Fluid package components of MSL 3.2.1"

      constant String substanceNames[:]= {""}
      "To express number and order of substances"         annotation (Evaluate=true);

      replaceable package Medium = Interfaces.SimpleChemicalMedium (substanceNames=substanceNames,
            redeclare package stateOfMatter = Interfaces.Incompressible)
      "Medium model"   annotation (choicesAllMatching=true);

      package StateOfMatter = Medium.stateOfMatter
      "State of matter of each chemical substance";
      constant StateOfMatter.SubstanceData substanceData[:] = Medium.substanceData
      "Definitions of all chemical substances";

      Modelica.Fluid.Interfaces.FluidPort_a fluid(redeclare package Medium = Medium)
      "Connector of Modelica.Fluid package"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      Interfaces.SubstancePorts_b substances[n]
      "All chemical substances of the solution"                                           annotation (Placement(transformation(
              extent={{-110,-40},{-90,40}}), iconTransformation(extent={{-110,-40},{
                -90,40}})));
      Interfaces.SolutionPort solution "Chemical solution"
        annotation (Placement(transformation(extent={{30,-40},{50,-20}}),
            iconTransformation(extent={{30,-40},{50,-20}})));
       Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
      "Heat port of the chemical solution"                                                               annotation (
          Placement(transformation(extent={{-50,-40},{-30,-20}}),iconTransformation(
              extent={{-50,-40},{-30,-20}})));

      Modelica.SIunits.MoleFraction x[n] "Mole fraction of the substance";

  protected
      constant Integer n=size(substanceNames,1) "Number of substances";
      Modelica.SIunits.MolarMass molarMass[n] "Molar mass of the substance";

      Modelica.SIunits.Temperature temperature(start=298.15)
      "Temperature of the solution";

      Modelica.SIunits.Pressure pressure(start=100000)
      "Pressure of the solution";

      Modelica.SIunits.ElectricPotential electricPotential(start=0)
      "Electric potential of the solution";

      Medium.ThermodynamicState actualStreamThermodynamicState
      "Thermodynamic state of solution inside stream";

      //  Modelica.SIunits.SpecificEnergy actualStreamSpecificEnthalpy
      //    "Specific Enthalpy of solution inside stream";

      Modelica.SIunits.SpecificEnergy actualStreamSpecificInternalEnergy
      "Specific Internal Energy of solution inside stream";

    equation
      //fluid
      fluid.p = pressure;
      fluid.m_flow * actualStream(fluid.Xi_outflow) + molarMass .* substances.q = zeros(n);
      fluid.Xi_outflow = x .* molarMass ./ (x*molarMass);

      //substances
      molarMass = StateOfMatter.molarMass(Medium.substanceData,temperature,pressure,electricPotential,solution.I);
      substances.u = StateOfMatter.chemicalPotentialPure(
        Medium.substanceData,
        temperature,
        pressure,
        electricPotential,
        solution.I)
        + Modelica.Constants.R*temperature*log(x .* StateOfMatter.activityCoefficient(Medium.substanceData,temperature,pressure,electricPotential,solution.I))
        + Modelica.Constants.F*electricPotential*StateOfMatter.chargeNumberOfIon(Medium.substanceData,temperature,pressure,electricPotential,solution.I);

      //energy balance
      //fluid.h_outflow = solution.H / solution.m;
      fluid.h_outflow = Medium.specificEnthalpy(Medium.setState_pTX(pressure, temperature, fluid.Xi_outflow, electricPotential, solution.I));

      actualStreamThermodynamicState = Medium.setState_phX(pressure, actualStream(fluid.h_outflow), actualStream(fluid.Xi_outflow), electricPotential, solution.I);

      //The first idea was to changne only the enthalpy as extensive energy, but it changes the temperature with change of volume during isobaric and isothermic conditions!!!
      //  actualStreamSpecificEnthalpy = Medium.specificEnthalpy(actualStreamThermodynamicState);
      //  heatPort.Q_flow = -actualStreamSpecificEnthalpy*fluid.m_flow;

      //As a result, whole internal energy must be changed during mass changes. This energy is balenced using heat port with solution.
      actualStreamSpecificInternalEnergy = Medium.specificInternalEnergy(actualStreamThermodynamicState);
      heatPort.Q_flow = -actualStreamSpecificInternalEnergy*fluid.m_flow;

      //solution aliasses
      temperature = solution.T;
      pressure = solution.p;
      electricPotential = solution.v;

      //do not affect solution at port:
      solution.dH = 0;
      solution.i = 0;
      solution.dV = 0;
      solution.Gj = 0;
      solution.nj = 0;
      solution.mj = 0;
      solution.Qj = 0;
      solution.Ij = 0;
      solution.Vj = 0;

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={Line(
              points={{-90,0},{90,0}},
              color={158,66,200},
              thickness=1)}));
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
              color={127,0,127})}),
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
              color={127,0,127})}),
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
              color={127,0,127})}),
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
              color={127,0,127})}),
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
              color={127,0,127})}),
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
              color={127,0,127})}),
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
      annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

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
        annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

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
        annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

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
        annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

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
        annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

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
          annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

          extends Interfaces.ConditionalKinetics(KC=1/(Modelica.Constants.R*298.15));

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
       Text(extent=  {{-160,110},{40,50}}, lineColor={172,72,218},   textString=  "%name")}),
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
       Text(extent=  {{-160,110},{40,50}}, lineColor={172,72,218},   textString=  "%name")}),
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

     replaceable function molarHeatCapacityCp
      "Molar heat capacity at constant pressure"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarHeatCapacity molarHeatCapacityCp
        "Molar heat capacity at constant pressure";
     end molarHeatCapacityCp;

      replaceable function molarHeatCapacityCv
      "Molar heat capacity at constant volume"
        extends Modelica.Icons.Function;
        input SubstanceData substanceData "Data record of substance";
        input Modelica.SIunits.Temperature T=298.15 "Temperature";
        input Modelica.SIunits.Pressure p=100000 "Pressure";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the substance";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        output Modelica.SIunits.MolarHeatCapacity molarHeatCapacityCv
        "Molar heat capacity at constant volume";
      end molarHeatCapacityCv;

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
        "Enthalpy of formation of the substance at 25 degC";

          parameter Modelica.SIunits.MolarEnergy DfG_25degC_1bar(displayUnit="kJ/mol")=0
        "Gibbs enerfy of formation of the substance at 25 degC,1bar";

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
           +(T-298.15)*(substanceData.Cp);
     end molarEnthalpyElectroneutral;

     redeclare function extends molarEntropyPure
      "Molar entropy of the pure substance"
     algorithm
       //molarEntropyPure := ((substanceData.DfH - substanceData.DfG_25degC_1bar)/298.15)
       //+ (substanceData.Cp+Modelica.Constants.R)*log(T/298.15);

         //Molar entropy:
         // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
         // - pressure shift: to reach the ideal gas equation at constant temperature Vm*dP = -T*dS (small amount of work)
         molarEntropyPure := (substanceData.Cp)*log(T/298.15) - Modelica.Constants.R*log(p/100000) + ((substanceData.DfH_25degC
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

     redeclare function extends molarHeatCapacityCp
      "Molar heat capacity of the substance at constant pressure"
     algorithm
         molarHeatCapacityCp := substanceData.Cp;
     end molarHeatCapacityCp;

     redeclare function extends molarHeatCapacityCv
      "Molar heat capacity of the substance at constant volume"
     algorithm
         molarHeatCapacityCv := substanceData.Cp - Modelica.Constants.R;
     end molarHeatCapacityCv;

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
        "Enthalpy of formation of the substance at 25 degC";

          parameter Modelica.SIunits.MolarEnergy DfG_25degC_1bar(displayUnit="kJ/mol")=0
        "Gibbs enerfy of formation of the substance at 25 degC,1bar";

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

     redeclare function extends molarHeatCapacityCp
      "Molar heat capacity of the substance at constant pressure"
     algorithm
         molarHeatCapacityCp := substanceData.Cp;
     end molarHeatCapacityCp;

     redeclare function extends molarHeatCapacityCv
      "Molar heat capacity of the substance at constant volume"
     algorithm
         molarHeatCapacityCv := substanceData.Cp;
     end molarHeatCapacityCv;

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
              lineColor={127,127,0},
              fillColor={127,127,0},
              fillPattern=FillPattern.Solid)}),
      Diagram(graphics={
       Text(extent={{-160,110},{40,50}},   lineColor={127,127,0},    textString=  "%name",
            fillColor={127,127,0},
            fillPattern=FillPattern.Solid),
                      Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={127,127,0},
              fillColor={127,127,0},
              fillPattern=FillPattern.Solid,
              lineThickness=1)}));
    end SolutionPort;

    partial model PartialSolution
    "Chemical solution as homogenous mixture of the substances (only pressure and electric potential are not defined)"

      Modelica.SIunits.Temperature T(start=298.15) "Temperature";

      Modelica.SIunits.Pressure p(start=100000) "Pressure";

      Modelica.SIunits.Volume volume(stateSelect=StateSelect.prefer)
      "Current volume of the solution";

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

    //   Modelica.SIunits.HeatFlowRate der_freeEnthalpy;
    initial equation

      freeInternalEnergy = 0;
      //freeEnthalpy + solution.Hj = 0;
      //der(freeEnthalpy) = -solution.dH;
    equation
      //internal energy
      der(freeInternalEnergy) = heatFromEnvironment + workFromEnvironment;

      heatFromEnvironment + workFromEnvironment = (-solution.dH) - solution.p*(-solution.dV) - volume*der(solution.p);
     // heatFromEnvironment + workFromEnvironment = der(freeEnthalpy) - solution.p*(-solution.dV) - volume*der(solution.p);
      //It is the same as: der(freeEnthalpy)=-solution.dH;

      //thermodinamics equations:
      freeInternalEnergy = freeEnthalpy - volume*p; // H=U+p*V
      freeGibbsEnergy = freeEnthalpy - T*freeEntropy; // G=H-T*S

      //aliases
      solution.p = p;
      solution.T = T;
      solution.G = freeGibbsEnergy;
    //  solution.H = freeEnthalpy;

      //  solution.U = freeInternalEnergy;
      solution.Q = charge;
      solution.V = volume;
    //  der_freeEnthalpy = der(freeEnthalpy);

      //Extensive properties of the solution:

      // The extensive quantities here have not the real physical flows.
      // They hack the Kirchhof's flow equation to be counted as the sum from all connected substances in the solution.

      //amount of substances
      solution.n + solution.nj = 0; //total amount of solution is the sum of amounts of each substance

      //mass of substances
      solution.m + solution.mj = 0; //total mass of solution is the sum masses of each substance

      //free Gibs energy
      solution.G + solution.Gj = 0; //total free Gibbs energy of solution is the sum of free Gibbs energies of each substance

      //free enthalpy
    //  solution.H + solution.Hj = 0;  //total free enthalpy of solution is the sum of enthalpies of each substance

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
</html>"));
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
      annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

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
      annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

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
      annotation(Evaluate=true, HideResult=true, choices(checkbox=true),Dialog(group="Conditional inputs"));

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

    package SimpleChemicalMedium
    "Mixture of chemical substances as homogenous solution in one state of matter"
      extends Modelica.Media.Interfaces.PartialMedium(
        final mediumName="ChemicalSolution",
        final singleState=false,
        final reducedX=false,
        final fixedX=false,
        Temperature(
          min=273,
          max=450,
          start=298));

      replaceable package stateOfMatter =
                              Interfaces.Incompressible  constrainedby
      Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
         annotation (choicesAllMatching = true);

      // Provide medium constants here
      constant stateOfMatter.SubstanceData substanceData[nS]
      "Definition of the substances"
         annotation (choicesAllMatching = true);

      // Simplified properties - default values are for water 25degC from Modelica.Media.CompressibleLiquids.LinearWater_pT_Ambient
      constant DynamicViscosity eta_const = 8.9e-4 "Constant dynamic viscosity";

      redeclare model extends BaseProperties(final standardOrderComponents=true)
      "Base properties of medium"

       input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the solution";
       input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

       Modelica.SIunits.MolarMass molarMass[nS] "Molar mass";
       Modelica.SIunits.MolarVolume molarVolume[nS] "Molar volume";
       Modelica.SIunits.MolarEnthalpy molarEnthalpy[nS] "Molar enthalpy";
       Modelica.SIunits.MolarHeatCapacity molarHeatCapacity[nS]
        "Molar heat capacity";

       SpecificHeatCapacity Cp "Specific heat capacity at constant pressure";
      equation

        molarMass = stateOfMatter.molarMass(substanceData,T,p,v,I);
        molarVolume = stateOfMatter.molarVolume(substanceData,T,p,v,I);
        molarEnthalpy = stateOfMatter.molarEnthalpy(substanceData,T,p,v,I);
        molarHeatCapacity = stateOfMatter.molarHeatCapacityCp(substanceData,T,p,v,I);

        MM = sum( molarMass ./ X);  // sum(MMj*nj)/sum(nj)
        d = 1 / sum( molarVolume .* X ./ molarMass);  //sum(MMj*nj)/sum(Vmj*nj), where xj=nj/nT, sum(Xj)=1
        Cp = sum( molarHeatCapacity .* X ./ molarMass); //sum(cpj*nj)/sum(MMj*nj)
        h = sum( molarEnthalpy .* X ./ molarMass); //sum(Hmj*nj)/sum(MMj*nj)
        u = h - p/d;
        R = Modelica.Constants.R/MM;
        state.p = p;
        state.T = T;
        state.X = X;
        state.v = v; //electric potential
        state.I = I; //ionic strength
      end BaseProperties;

      redeclare record ThermodynamicState
      "A selection of variables that uniquely defines the thermodynamic state"
        extends Modelica.Icons.Record;
        AbsolutePressure p "Absolute pressure of medium";
        Temperature T "Temperature of medium";

        Modelica.SIunits.MassFraction X[nS] "Mass fractions";

        Modelica.SIunits.ElectricPotential v
        "Electric potential of the solution";
        Modelica.SIunits.MoleFraction I "Ionic strengh (mole fraction based)";
        annotation (Documentation(info="<html>

</html>"));
      end ThermodynamicState;

      redeclare function setState_pTX
      "Return thermodynamic state from p, T, and X or Xi"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input Temperature T "Temperature";
      input MassFraction X[:]=reference_X "Mass fractions";
      input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the solution";
      input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        //  input stateOfMatter.SubstanceData substanceData[nS] = substanceData       "Constant data of substances";
      output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := ThermodynamicState(p=p, T=T, X=X, v=v, I=I);
        annotation (Documentation(info="<html>

</html>"));
      end setState_pTX;

      redeclare function setState_phX
      "Return thermodynamic state from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[nS]=reference_X "Mass fractions";
        input Modelica.SIunits.ElectricPotential v=0
        "Electric potential of the solution";
        input Modelica.SIunits.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
        //    input stateOfMatter.SubstanceData substanceData[nS] = substanceData     "Constant data of substances";
        output ThermodynamicState state "Thermodynamic state record";
    protected
       Modelica.SIunits.MolarMass molarMass[nS] "Molar mass";
       Modelica.SIunits.MolarEnthalpy molarEnthalpyT0[nS]
        "Molar enthalpy at temperature T0";
       Modelica.SIunits.MolarHeatCapacity molarHeatCapacity[nS]
        "Molar heat capacity";

       SpecificHeatCapacity Cp "Specific heat capacity at constant pressure";
       Modelica.SIunits.SpecificEnthalpy hT0
        "Molar enthalpy of solution at temperature T0";
       constant Temperature T0=298;
      algorithm
        molarMass := stateOfMatter.molarMass(substanceData,T0,p,v,I);
        molarEnthalpyT0 :=stateOfMatter.molarEnthalpy(
            substanceData,
            T0,
            p,
            v,
            I);
        hT0 := sum(molarEnthalpyT0 .* X ./ molarMass);
        molarHeatCapacity := stateOfMatter.molarHeatCapacityCp(substanceData,T0,p,v,I);
        Cp := sum( molarHeatCapacity .* X ./ molarMass);
        state := ThermodynamicState(p=p, T=T0 + (h-hT0)/Cp, X=X, v=v, I=I);
      end setState_phX;

      //TODO: make model of dynamic viscosity from substance properties!!!
      redeclare function extends dynamicViscosity
      "Return constant dynamic viscosity"
      algorithm
        eta := eta_const;
      end dynamicViscosity;

      redeclare function extends pressure "Return pressure"

      algorithm
        p := state.p;
      end pressure;

      redeclare function extends temperature "Return temperature"

      algorithm
        T := state.T;
      end temperature;

      redeclare function extends density "Return density"
        //  input stateOfMatter.SubstanceData substanceData[nS] = substanceData       "Constant data of substances";
    protected
       Modelica.SIunits.MolarMass molarMass[nS] "Molar mass";
       Modelica.SIunits.MolarVolume molarVolume[nS] "Molar volume";
      algorithm
        molarMass := stateOfMatter.molarMass(substanceData,state.T,state.p,state.v,state.I);
        molarVolume := stateOfMatter.molarVolume(substanceData,state.T,state.p,state.v,state.I);

        d := 1 / sum( molarVolume .* state.X ./ molarMass); //sum(MMj*nj)/sum(Vmj*nj), where xj=nj/nT, sum(Xj)=1

      end density;

      redeclare function extends specificEnthalpy "Return specific enthalpy"
        //    input stateOfMatter.SubstanceData substanceData[nS] = substanceData       "Constant data of substances";
      algorithm
        h := sum( stateOfMatter.molarEnthalpy(substanceData,state.T,state.p,state.v,state.I) .* state.X ./ stateOfMatter.molarMass(substanceData,state.T,state.p,state.v,state.I));
      end specificEnthalpy;

      redeclare function extends specificEntropy
      "Return the specific entropy from the thermodynamic state"
        //    input stateOfMatter.SubstanceData substanceData[nS] = substanceData     "Constant data of substances";
    protected
       Modelica.SIunits.MolarMass molarMass[nS] "Molar mass";
       Modelica.SIunits.MolarEntropy molarEntropy[nS] "Molar entropies";
       Modelica.SIunits.MoleFraction x[nS] "Mole fractions";
       Modelica.SIunits.ChemicalPotential electrochemicalPotential[nS]
        "electrochemical potentials";
      algorithm
        molarMass := stateOfMatter.molarMass(substanceData,state.T,state.p,state.v,state.I);

        x := (state.X ./ molarMass) / sum(state.X ./ molarMass);
        electrochemicalPotential :=
          stateOfMatter.electroChemicalPotentialPure(substanceData,state.T,state.p,state.v,state.I) +
          Modelica.Constants.R * state.T * log(x.*stateOfMatter.activityCoefficient(substanceData,state.T,state.p,state.v,state.I))
          + Modelica.Constants.F * state.v * stateOfMatter.chargeNumberOfIon(substanceData,state.T,state.p,state.v,state.I);

        molarEntropy := stateOfMatter.molarEntropy(electrochemicalPotential,substanceData,state.T,state.p,state.v,state.I);

        s := sum( molarEntropy .* state.X ./ molarMass);
      end specificEntropy;

      redeclare function extends specificInternalEnergy
      "Return the specific internal energy from the thermodynamic state"
        //    input stateOfMatter.SubstanceData substanceData[nS] = substanceData     "Constant data of substances";
      algorithm
        u := specificEnthalpy(state) - state.p/density(state);
      end specificInternalEnergy;

      redeclare function extends specificGibbsEnergy
      "Return specific Gibbs energy from the thermodynamic state"
        extends Modelica.Icons.Function;
        //    input stateOfMatter.SubstanceData substanceData[nS] = substanceData
      //    "Constant data of substances";
      algorithm
        g := specificEnthalpy(state) - state.T*specificEntropy(state);
      end specificGibbsEnergy;

      redeclare function extends specificHelmholtzEnergy
      "Return specific Helmholtz energy from the thermodynamic state"
        extends Modelica.Icons.Function;
      algorithm
        f := specificInternalEnergy(state) - state.T*specificEntropy(state);
      end specificHelmholtzEnergy;

      redeclare function extends specificHeatCapacityCp
      "Return specific heat capacity at constant pressure"
        //    input stateOfMatter.SubstanceData substanceData[nS] = substanceData
        //     "Constant data of substances";
      algorithm
        cp := sum( stateOfMatter.molarHeatCapacityCp(substanceData,state.T,state.p,state.v,state.I) .* state.X ./ stateOfMatter.molarMass(substanceData,state.T,state.p,state.v,state.I));

        annotation (Documentation(info="<html>

</html>"));
      end specificHeatCapacityCp;

      redeclare function extends specificHeatCapacityCv
      "Return specific heat capacity at constant volume"
       // input stateOfMatter.SubstanceData substanceData[nS] = substanceData
         // "Constant data of substances";
      algorithm
        cv := sum( stateOfMatter.molarHeatCapacityCv(substanceData,state.T,state.p,state.v,state.I) .* state.X ./ stateOfMatter.molarMass(substanceData,state.T,state.p,state.v,state.I));

        annotation (Documentation(info="<html>

</html>"));
      end specificHeatCapacityCv;

      annotation (Documentation(info="<HTML>
<p>
This package is a <b>template</b> for <b>new medium</b> models. For a new
medium model just make a copy of this package, remove the
\"partial\" keyword from the package and provide
the information that is requested in the comments of the
Modelica source.
</p>
</HTML>"));
    end SimpleChemicalMedium;
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
versionDate="2015-09-15",
dateModified = "2015-09-15 17:14:41Z",
conversion(
  from(version="1.1.0alpha", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.0_to_1.1.mos"),
  from(version="1.0.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.0_to_1.1.mos")),
uses(Modelica(version="3.2.1")),
  Documentation(revisions="<html>
<p>Copyright (c) 2008-2015, Marek Matej&aacute;k, Charles University in Prague </p>
<p>All rights reserved. </p>
<p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: </p>
<ol>
<li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. </li>
<li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. </li>
<li>Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. </li>
</ol>
<p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS &QUOT;AS IS&QUOT; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
</html>", info="<html>
<p>During each electro-chemical process an <a href=\"modelica://Chemical.Components.Substance\">electro-chemical potential</a> of the substances is equilibrating and all thermodynamical properties of the homogenous chemical solutions are evaluated. </p>
<p>Processes: chemical reactions, gas dissolution, diffusion, membrane transports, osmotic fluxes, electrochemical cells, electrodes, ..</p>
<p>Please see the <a href=\"modelica://Chemical.UsersGuide.Overview\">overview</a>.</p>
</html>"));
end Chemical;
