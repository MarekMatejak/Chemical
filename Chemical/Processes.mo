within Chemical;
package Processes "Undirected process package"
  model Reaction "Chemical Reaction"
    extends Chemical.Processes.Internal.PartialReactionWithProductsDefinition;
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = if nP==0 then SolutionChoice.Parameter else SolutionChoice.FirstSubstrate);

    import Chemical.Utilities.Types.SolutionChoice;
    import Chemical.Processes.Internal.Kinetics;

    replaceable function uDiff = Kinetics.traditionalPotentialDiff
      constrainedby
        Kinetics.partialPotentialDiff "Electro-chemical potential difference function"
      annotation(choicesAllMatching=true, Dialog(tab="Advanced"), Documentation(info="<html>
    <p>Electro-chemical potential loss function used in the diffusion.</p>
    </html>"));
    extends Chemical.Interfaces.ConditionalKinetics(k_forward=1);


     Chemical.Interfaces.Definition processDefinition;

  equation

    processDefinition = p*products.definition - s*substrates.definition;

    //chemical kinetics

    du = uDiff(rr,kf*Sx_fore,(kf/Kx)*Px_rear,solutionState,n_flow_coef_reg);


    //chemical solution and its propagation
    if nS>0 then
      connect(substrates[1].solution_forwards,inputSubstrateSolution);
      products.solution_forwards = fill(solutionState,nP);
      substrates.solution_rearwards = fill(solutionState,nS);
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
            extent={{-146,-72},{142,-42}},
            textColor={128,0,255},
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
</html>",   info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichiometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
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
<td><p>stoichiometric coefficients of i-th substance</p></td>
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

  model Diffusion "Solute diffusion"
    extends Icons.Diffusion;
    extends Chemical.Interfaces.SISO;
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = SolutionChoice.Parameter);
    extends Chemical.Interfaces.ConditionalKinetics
                                          (k_forward=1);

    import Chemical.Utilities.Types.SolutionChoice;



    replaceable function uDiff =
        Chemical.Processes.Internal.Kinetics.traditionalPotentialDiff
      constrainedby
        Internal.Kinetics.partialPotentialDiff "Electro-chemical potential difference function"
      annotation(choicesAllMatching=true, Dialog(tab="Advanced"), Documentation(info="<html>
    <p>Electro-chemical potential loss function used in the diffusion.</p>
    </html>"));

  equation
    fore.definition = rear.definition;
    fore.solution_forwards = solutionState;
    rear.solution_rearwards = solutionState;

    connect(rear.solution_forwards,inputSubstrateSolution);


    du = uDiff(n_flow,kf*Sx_fore,(kf/Kx)*Px_rear,solutionState,n_flow_coef_reg);


     annotation ( Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
  end Diffusion;

  model GasSolubility "Henry's law of gas solubility into liquid."
    extends Icons.GasSolubility;

    extends Chemical.Interfaces.SISO;
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = SolutionChoice.Parameter);
    extends Chemical.Interfaces.ConditionalKinetics(k_forward=1);

    import Chemical.Utilities.Types.SolutionChoice;
    import Chemical.Utilities.Types.FirstProductChoice;

    parameter FirstProductChoice productFrom=Chemical.Utilities.Types.FirstProductChoice.Substance "Choice of products definition"
      annotation (
    //Chemical.Utilities.Types.FirstProductChoice.Process "Choice of products definition"
                  HideResult=true, Dialog(group="Product definition"));

    parameter Chemical.Interfaces.Definition product = Chemical.Substances.Liquid.Unknown "Product definitions"
      annotation (choicesAllMatching=true, Dialog(group="Product definition", enable=(productFrom ==Chemical.Utilities.Types.FirstProductChoice.fromParameter)));

     parameter Chemical.Interfaces.Definition process = Chemical.Interfaces.processData(1)
     "Process definition"
        annotation (Dialog(enable=(productFrom==Chemical.Utilities.Types.FirstProductChoice.fromProcessEnergies)));

    replaceable function uDiff =
        Chemical.Processes.Internal.Kinetics.traditionalPotentialDiff
      constrainedby
        Internal.Kinetics.partialPotentialDiff "Electro-chemical potential difference function"
      annotation(choicesAllMatching=true, Dialog(tab="Advanced"), Documentation(info="<html>
    <p>Electro-chemical potential loss function used in the diffusion.</p>
    </html>"));


  equation

    if (productFrom ==FirstProductChoice.Substance)  then
      fore.definition = product;
    else
      fore.definition = rear.definition + process;
    end if;

    fore.solution_forwards = solutionState;
    rear.solution_rearwards = solutionState;
    connect(rear.solution_forwards,inputSubstrateSolution);

    du = uDiff(n_flow,kf*Sx_fore,(kf/Kx)*Px_rear,solutionState,n_flow_coef_reg);


    annotation (
     Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Gaseous substance dissolution in liquid (Henry&apos;s law, Raoult&apos;s law, Nernst dissolution in one). </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K<sub>H</sub> =x<sub>L</sub> / x<sub>g</sub>&nbsp;</p></td>
<td><p>Henry&apos;s coefficient, Raoult&apos;s coefficient</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>G = &Delta;<sub>f</sub>G<sub>L </sub>- &Delta;<sub>f</sub>G<sub>g </sub>= &Delta;<sub>sol</sub>H - T&middot;&Delta;<sub>sol</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K<sub>H</sub>&middot; (f<sub>L</sub> / f<sub>g</sub>)) </p></td>
<td><p>molar Gibb&apos;s energy of the dissolution</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>H = &Delta;<sub>f</sub>H<sub>L </sub>- &Delta;<sub>f</sub>H<sub>g</sub></p></td>
<td><p>molar enthalpy of the dissolution</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>S = &Delta;<sub>f</sub>S<sub>L</sub> - &Delta;<sub>f</sub>S<sub>g</sub> = <a href=\"modelica://Modelica.Constants\">k</a>&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(&Delta;<sub>sol</sub>&omega;) </p></td>
<td><p>molar entropy of the dissolution</p></td>
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

  model Membrane "Passive transport of the substance through semipermeable membrane"
    extends Icons.Membrane;

    extends Chemical.Interfaces.SISO;
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = SolutionChoice.Parameter);
    extends Chemical.Interfaces.ConditionalKinetics
                                          (k_forward=1);

    import Chemical.Utilities.Types.SolutionChoice;



    replaceable function uDiff =
        Chemical.Processes.Internal.Kinetics.traditionalPotentialDiff
      constrainedby
        Internal.Kinetics.partialPotentialDiff "Electro-chemical potential difference function"
      annotation(choicesAllMatching=true, Dialog(tab="Advanced"), Documentation(info="<html>
    <p>Electro-chemical potential loss function used in the diffusion.</p>
    </html>"));

  equation
    fore.definition = rear.definition;
    fore.solution_forwards = solutionState;
    rear.solution_rearwards = solutionState;

    connect(rear.solution_forwards,inputSubstrateSolution);

    du = uDiff(n_flow,kf*Sx_fore,(kf/Kx)*Px_rear,solutionState,n_flow_coef_reg);


    annotation ( Documentation(info="<html>
<p><u><b><font style=\"color: #008000; \">Filtration through semipermeable membrane.</font></b></u></p>
<p>The penetrating particles are driven by electric and chemical gradient to reach Donnan&apos;s equilibrium.</p>
<p>If zero-flow Donnan&apos;s equilibrium is reached. </p>
</html>",
        revisions="<html>
<p><i>2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
         Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics={
          Text(
            extent={{-97,-12},{97,12}},
            textString="%name",
            textColor={128,0,255},
          origin={69,2},
          rotation=90)}));
  end Membrane;

  model Pump "Prescribed substance molar flow"
    extends Chemical.Interfaces.SISO;
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = SolutionChoice.Parameter);
    extends Interfaces.ConditionalSubstanceFlow;

    import Chemical.Utilities.Types.SolutionChoice;


  equation
    fore.definition = rear.definition;
    fore.solution_forwards = solutionState;
    rear.solution_rearwards = solutionState;

    connect(rear.solution_forwards,inputSubstrateSolution);




    n_flow = q;

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
            textColor={128,0,255},
            origin={0,-72},
            rotation=360,
            textString="%name")}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end Pump;

  model ForwardReaction "Chemical Reaction"
    extends Chemical.Processes.Internal.PartialReactionWithProductsDefinition(final process, final firstProductFrom=Chemical.Utilities.Types.FirstProductChoice.Substance);
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = if nS==0 then SolutionChoice.Parameter else SolutionChoice.FirstSubstrate);
    extends Chemical.Interfaces.ConditionalKinetics(k_forward=1);


    import Chemical.Utilities.Types.SolutionChoice;


  equation

    rr = kf*Sx_fore;

    /*if (nS>0) then
    du_fore = -du_rear;
  end if;*/

    //chemical solution and its propagation
    connect(substrates[1].solution_forwards,inputSubstrateSolution);
    if nS>0 then
      products.solution_forwards = fill(solutionState,nP);
      substrates.solution_rearwards = fill(solutionState,nS);
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
            extent={{-160,-70},{160,-40}},
            textColor={128,0,255},
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
</html>",   info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichiometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
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
<td><p>stoichiometric coefficients of i-th substance</p></td>
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

  model Stream "Flow of whole solution"
    extends Chemical.Interfaces.SISO;
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = SolutionChoice.Parameter);
    extends Boundaries.Internal.ConditionalSolutionFlow;

    import Chemical.Utilities.Types.SolutionChoice;


  equation
    fore.definition = rear.definition;
    fore.solution_forwards = solutionState;
    rear.solution_rearwards = solutionState;

    connect(rear.solution_forwards,inputSubstrateSolution);

    n_flow = Sx_fore * (solutionState.n/solutionState.V) * volumeFlow;


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
            textColor={128,0,255},
            origin={2,-74},
            rotation=180)}),
      Documentation(revisions="<html>
<p><i>2009-2018 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
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

  model Allostery "Change of macromolecule state - providing specific subunit forms (e.g. relaxed<->tensed deoxygenated hemoglobin change)"
    extends Chemical.Processes.Internal.PartialSpeciationWithSubunitsDefinition;
   // extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = if nP==0 then SolutionChoice.Parameter else SolutionChoice.FirstSubstrate);

    //import Chemical.Utilities.Types.SolutionChoice;
    import Chemical.Processes.Internal.Kinetics;

    replaceable function uDiff = Kinetics.traditionalPotentialDiff
      constrainedby
        Kinetics.partialPotentialDiff "Electro-chemical potential difference function"
      annotation(choicesAllMatching=true, Documentation(info="<html>
    <p>Electro-chemical potential loss function used in the diffusion.</p>
    </html>"));
    extends Chemical.Interfaces.ConditionalKinetics(k_forward=1);



  equation

    //chemical kinetics
    du = uDiff(rr,kf*Sx_fore,(kf/Kx)*Px_rear,solutionState,n_flow_coef_reg);


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
            textColor={128,0,255},
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
</html>",   info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichiometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
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
<td><p>stoichiometric coefficients of i-th substance</p></td>
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
  end Allostery;

  model MichaelisMenten "Enzyme catalysed reaction using Michaelis-Menten kinetics"
    extends Chemical.Icons.EnzymeKinetics;
    extends Chemical.Processes.Internal.PartialReactionWithProductsDefinition(final process, final firstProductFrom=Chemical.Utilities.Types.FirstProductChoice.Substance);
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = if nS==0 then SolutionChoice.Parameter else SolutionChoice.FirstSubstrate);


    import Chemical.Utilities.Types.SolutionChoice;

    parameter Real k_cat;
    parameter Real Km;

    Modelica.Blocks.Interfaces.RealInput e0 "Initial enzyme concentration" annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-42,100}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-42,100})));
  equation


    rr = k_cat*e0*Sx_fore/(Km+Sx_fore);

    //chemical solution and its propagation
    connect(substrates[1].solution_forwards,inputSubstrateSolution);
    if nS>0 then
      products.solution_forwards = fill(solutionState,nP);
      substrates.solution_rearwards = fill(solutionState,nS);
    end if;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-160,-70},{160,-40}},
            textColor={128,0,255},
          textString="%name")}),
      Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichiometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
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
<td><p>stoichiometric coefficients of i-th substance</p></td>
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
  end MichaelisMenten;

  model Antagonist
    parameter Real kf,kb,n_complex_start;

    Reaction binding(
      process=Chemical.Interfaces.processData(kf/kb),
      k_forward=kf,
      nP=1,
      nS=2) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={22,-12})));
    Boundaries.Substance substance(
      useRear=true,
      preferMass=false,
      amountOfSubstance_start(displayUnit="umol") = n_complex_start) annotation (Placement(transformation(extent={{56,-22},{76,-2}})));
    Topology.JunctionRFF2 junctionRFF2 annotation (Placement(transformation(extent={{-28,70},{-8,50}})));
    Topology.JunctionRFF2 junctionRFF1 annotation (Placement(transformation(extent={{-26,-90},{-6,-70}})));
    Interfaces.Fore fore annotation (Placement(transformation(extent={{92,50},{112,70}})));
    Interfaces.Rear rear annotation (Placement(transformation(extent={{-110,50},{-90,70}})));
    Interfaces.Fore foreLigand annotation (Placement(transformation(extent={{90,-90},{110,-70}})));
    Interfaces.Rear rearLigand annotation (Placement(transformation(extent={{-110,-90},{-90,-70}})));
  equation
    connect(binding.products[1], substance.rear) annotation (Line(
        points={{32,-12},{56,-12}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionRFF2.foreA, binding.substrates[1])
      annotation (Line(
        points={{-18,50},{-18,-12.25},{12,-12.25}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionRFF1.foreA, binding.substrates[2])
      annotation (Line(
        points={{-16,-70},{-16,-11.75},{12,-11.75}},
        color={158,66,200},
        thickness=0.5));
    connect(rear, junctionRFF2.rear) annotation (Line(
        points={{-100,60},{-28,60}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionRFF2.foreB, fore) annotation (Line(
        points={{-8,60},{102,60}},
        color={158,66,200},
        thickness=0.5));
    connect(rearLigand, junctionRFF1.rear) annotation (Line(
        points={{-100,-80},{-26,-80}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionRFF1.foreB, foreLigand) annotation (Line(
        points={{-6,-80},{100,-80}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,60},{104,60}},
            color={158,66,200},
            thickness=1),
          Line(
            points={{-100,-80},{100,-80}},
            color={158,66,200},
            thickness=1),
          Rectangle(
            extent={{-40,40},{40,30}},
            lineColor={0,0,0},
            lineThickness=1,
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-6,34},{4,-60}},
            lineColor={0,0,0},
            lineThickness=1,
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),                      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end Antagonist;

  package Internal "Internals package for Processes"
    extends Modelica.Icons.InternalPackage;

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
    partial model PartialReactionWithProductsDefinition "Chemical Reaction"
      extends Chemical.Processes.Internal.PartialReaction;

      import Chemical.Substances.Liquid;
      import Chemical.Substances.Gas;
      import Chemical.Substances.Solid;
      import Chemical.Substances.Aqueous;

      import FirstProductChoice =
             Chemical.Utilities.Types.FirstProductChoice;

      parameter FirstProductChoice firstProductFrom=FirstProductChoice.Process  "First product definition comes from?"
          annotation(HideResult=true, Dialog(group="Products definitions"));

      parameter Chemical.Interfaces.Definition firstProduct=dropOfCommons.DefaultSubstance "First product definition as Substance"
        annotation (choicesAllMatching=true, Dialog(enable=(firstProductFrom == FirstProductChoice.Substance), group="Products definitions"));


      // for unknown reason this array can not be empty in Dymola 2025 - stupid fix: max(1,nP-1) instead of max(0,nP-1) is used as size of the array
      parameter Chemical.Interfaces.Definition nextProducts[:]=fill(((s*ones(nS))/(p*ones(nP)))*dropOfCommons.DefaultSubstance, max(1,nP-1)) "Definitions of next products"
        annotation (choicesAllMatching=true, Dialog(enable=(nP > 1), group="Products definitions"));


      parameter Chemical.Interfaces.Definition process = Chemical.Interfaces.processData(1)
       "Process definition"
       annotation (Dialog(enable=(firstProductFrom == FirstProductChoice.Process)));

    /*
  parameter Real K=1 "Process dissociation constant at 25°C,1bar"
   annotation (HideResult=true, Dialog(group="Process definition",enable=(firstProductFrom == FirstProductChoice.ProcessProperties)));

  parameter Modelica.Units.SI.MolarEnergy dH=0 "Process molar enthalpy change at 25°C,1bar"
   annotation (HideResult=true,Dialog(group="Process definition",enable=(firstProductFrom == FirstProductChoice.ProcessProperties)));

  parameter Modelica.Units.SI.MolarHeatCapacity dCp=0 "Process molar heat capacity change at 25°C,1bar"
   annotation (HideResult=true,Dialog(group="Process definition",enable=(firstProductFrom == FirstProductChoice.ProcessProperties)));

  parameter Modelica.Units.SI.SpecificVolume dVs=0 "Process specific volume change at 25°C,1bar [L/g]"
   annotation (HideResult=true,Dialog(group="Process definition",enable=(firstProductFrom == FirstProductChoice.ProcessProperties)));
*/
    equation

      if (nP>1) then
        for i in 1:nP-1 loop
           products[i+1].definition = nextProducts[i];
        end for;
      end if;

      if (nP>0)  then

        if (firstProductFrom == FirstProductChoice.Substance)  then
          products[1].definition = firstProduct;
        elseif (nP>1) then
          for i in 1:1 loop //this stupid loop is only for disabling check error about empty array p[:] in Dymola 2025
            products[1].definition =
              (1/p[i]) * (s*substrates.definition + process - p[2:end]*nextProducts);
          end for;

        else
          for i in 1:1 loop //this stupid loop is only for disabling check error about empty array p[:] in Dymola 2025
            products[1].definition = (1/p[i])*(s*substrates.definition + process);
          end for;
        end if;
      end if;



      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}})),
        Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichiometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
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
<td><p>stoichiometric coefficients of i-th substance</p></td>
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
    end PartialReactionWithProductsDefinition;

    partial model PartialReaction "Chemical Reaction"
      import Chemical;
      import Chemical.Utilities.Types.InitializationMethods;


      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
        annotation(HideResult=true, Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_coef_reg=dropOfCommons.n_flow_coef_reg "Regulation forward/backward flow tolerance of chemical kinectics - smallest significant rate near forward is (1-coef_reg)*q_f (resp. (1-coef_reg)*q_b for backward)";
      parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
        annotation(HideResult=true, Dialog(tab="Advanced"));
      parameter InitializationMethods initN_flow =Chemical.Utilities.Types.InitializationMethods.none  "Initialization method for n_flow"
        annotation(HideResult=true, Dialog(tab= "Initialization", group="Molar flow"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_0 = 0 "Initial value for n_flow"
        annotation(HideResult=true, Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.state)));
      parameter Utilities.Units.MolarFlowAcceleration n_acceleration_0 = 0 "Initial value for der(n_flow)"
        annotation(HideResult=true, Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.derivative)));
      parameter Modelica.Units.SI.Time TC = dropOfCommons.TC "Time constant for electro-chemical potential adaption" annotation (HideResult=true, Dialog(tab="Advanced"));
      parameter Utilities.Units.Inertance L = dropOfCommons.L "Inertance of the flow"
        annotation(HideResult=true, Dialog(tab="Advanced"));
      parameter Integer nS=0 "Number of substrate types"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));
      parameter Modelica.Units.SI.StoichiometricNumber s[nS]=ones(nS)
        "Stoichiometric coefficients for substrates"
        annotation (HideResult=true);
      parameter Integer nP=0 "Number of product types"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));
      parameter Modelica.Units.SI.StoichiometricNumber p[nP]=ones(nP)
        "Stoichiometric numbers for products"
        annotation (HideResult=true);
      Modelica.Units.SI.MolarFlowRate rr(stateSelect=n_flowStateSelect) "Reaction molar flow rate";
      Chemical.Interfaces.Rear substrates[nS] annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-100,0}), iconTransformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-100,0})));
      Chemical.Interfaces.Fore products[nP] annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={100,0}), iconTransformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={100,0})));
     protected
      Modelica.Units.SI.MolarEnthalpy h_fore_mix, h_rear_mix;
      Real Sx_fore,Px_rear,Kx;
      Modelica.Units.SI.ChemicalPotential uPure_substrates[nS];
      Modelica.Units.SI.ChemicalPotential uPure_products[nP];
      outer DropOfCommons dropOfCommons;
      Modelica.Units.SI.ChemicalPotential substrates_r_intern[nS]=Chemical.Utilities.Internal.regStep(
                substrates.n_flow,
                s.*(p*products.r)./(s*ones(nS)),
                0,
                n_flow_reg);
      Modelica.Units.SI.ChemicalPotential products_r_intern[nP]=Chemical.Utilities.Internal.regStep(
                products.n_flow,
                p.*(s*substrates.r)./(p*ones(nP)),
                0,
                n_flow_reg);
      Real _rR[nS+nP];
      Real _uIn[nS+nP];
      Real _uOut[nS+nP];
      Real _qIn[nS+nP];
      Real _uS[nS],_uP[nP];
      Real du;
    initial equation
      if initN_flow == InitializationMethods.state then
        rr = n_flow_0;
      elseif initN_flow == InitializationMethods.derivative then
        der(rr) = n_acceleration_0;
      elseif initN_flow == InitializationMethods.steadyState then
        der(rr) = 0;
      end if;
    equation
      rr*s = substrates.n_flow;
      rr*p = -products.n_flow;
      s*(der(rr)*L) =  substrates.r - _rR[1:nS];
      -p*(der(rr)*L) =  products.r - _rR[nS+1:end];
      for i in 1:nS loop
        _uIn[i] = substrates[i].state_forwards.u;
        _uOut[i] = (_uIn*_qIn - _uIn[i]*_qIn[i])/((_qIn*ones(nS+nP)-_qIn[i])*(s*ones(nS)));
        _qIn[i] = max(substrates[i].n_flow,n_flow_reg);
        Chemical.Utilities.Internal.regStep(substrates[i].n_flow,_uIn[i],_uOut[i],n_flow_reg) + _rR[i] = _uS[i];
        substrates[i].state_rearwards.u = _uOut[i];
        uPure_substrates[i] = Chemical.Interfaces.Properties.electroChemicalPotentialPure(substrates[i].definition, substrates[i].solution_forwards);
      end for;
      for i in 1:nP loop
        _uIn[i+nS] = products[i].state_rearwards.u;
        _uOut[i+nS] = (_uIn*_qIn - _uIn[i+nS]*_qIn[i+nS])/((_qIn*ones(nS+nP)-_qIn[i+nS])*(p*ones(nP)));
        _qIn[i+nS] = max(products[i].n_flow,n_flow_reg);
        Chemical.Utilities.Internal.regStep(products[i].n_flow,_uIn[i+nS],_uOut[i+nS],n_flow_reg) + _rR[i+nS] = _uP[i];
        products[i].state_forwards.u = _uOut[i+nS];
        uPure_products[i] = Chemical.Interfaces.Properties.electroChemicalPotentialPure( products[i].definition,  products[i].solution_rearwards);
      end for;
      du = p*_uP - s*_uS;
      Px_rear = exp((p * ((products.state_rearwards.u - uPure_products)./(Modelica.Constants.R*products.solution_rearwards.T))));
      Sx_fore = exp(s * ((substrates.state_forwards.u - uPure_substrates)./(Modelica.Constants.R*substrates.solution_forwards.T)));
      Kx = exp( ((s * ((uPure_substrates)./(Modelica.Constants.R*substrates.solution_forwards.T))) - (p * ((uPure_products)./(Modelica.Constants.R*products.solution_rearwards.T)))));
      products.state_forwards.h = h_fore_mix*ones(nP);
      substrates.state_rearwards.h = h_rear_mix*ones(nS);
      if noEvent((nS>0) and (abs(substrates.n_flow*ones(nS))>Modelica.Constants.eps))  then
        h_rear_mix*(substrates.n_flow*ones(nS)) + products.n_flow*products.state_rearwards.h = 0;
      else
        h_rear_mix = 0;
      end if;
      if noEvent((nP>0) and (abs(products.n_flow*ones(nP))>Modelica.Constants.eps)) then
        h_fore_mix*(products.n_flow*ones(nP)) + substrates.n_flow*substrates.state_forwards.h = 0;
      else
        h_fore_mix = 0;
      end if;
        annotation(HideResult=true, Dialog(tab="Advanced"),
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}})),
        Documentation(revisions="<html>
<p><i>2013-2025 by </i>Marek Mateják </p>
</html>",     info="<html>
<h4><span style=\"color: #008000\">Notations</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>A<sub>i</sub></p></td>
<td><p>i-th substance</p></td>
</tr>
<tr>
<td><p>v<sub>i</sub></p></td>
<td><p>stoichiometric coefficients of i-th substance</p></td>
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

    package Kinetics

      partial function partialPotentialDiff "Interface for potential difference functions"
        extends Modelica.Icons.Function;

        input Modelica.Units.SI.MolarFlowRate q "Molar flow rate of the process";
        input Modelica.Units.SI.MolarFlowRate q_f "Forward rate (e.g. k_forward*x_substrates)";
        input Modelica.Units.SI.MolarFlowRate q_b "Backward rate (e.g. k_backward*x_products)";
        input Chemical.Interfaces.SolutionState solutionState "Solution state";
        input Real coef_reg = 1e-5 "Forward/backward flow tolerance - smallest significant rate near forward is (1-coef_reg)*q_f (resp. (1-coef_reg)*q_b for backward)";

        output Modelica.Units.SI.ChemicalPotential potentialDiff "Electro-chemical potential difference in chemical process";

        annotation(Inline=true, smoothOrder=100,
          Documentation(info="<html>
    <p>Interface definition for a potential loss in a chemical process. Inputs are information about flow condition and the chemical solution state, output is the electro-chemical potential drop.</p>
</html>"));
      end partialPotentialDiff;

      function traditionalPotentialDiff "Traditional potential difference function trough forward/backward rate"
        extends partialPotentialDiff;

      algorithm
        potentialDiff :=
          if (q > (1-coef_reg)*q_f) then Modelica.Constants.R*solutionState.T*log(coef_reg)
          elseif (q > 0) then Modelica.Constants.R*solutionState.T*log(abs(1-q/q_f))
          elseif (q < -(1-coef_reg)*q_b)  then -Modelica.Constants.R*solutionState.T*log(coef_reg)
          else -Modelica.Constants.R*solutionState.T*log(abs(1+q/q_b));
          annotation(Dialog(enable=true),
                    Documentation(info="<html>
<p>
This Gibbs energy loss du to reach chemical process molar flow q: 
</p>
<blockquote><pre>
 = kf * xA - kb * xB
</pre></blockquote>
<p>
where kf is the forward rate coefficient and kb is backward rate coefficient of chemical process;
xA is mole fraction os substrate and xB is mole fraction of product.
And K = kf/kb = xB/xA is a dissociation coefficient of the chemical process.
</p>
</html>"));

      end traditionalPotentialDiff;

      function linearPotentialDiff "Linear potential difference function"
        extends partialPotentialDiff;

        input Real kC(unit="mol2/(J.s)") = 1 "Linear factor"
          annotation(Dialog(enable=true));


      algorithm
        potentialDiff := -q/kC;

        annotation (Documentation(info="<html>
<p>
This Gibbs energy loss is linear in the molarflow with the linear factor kC: 
</p>
<blockquote><pre>
du := n_flow/kC;
</pre></blockquote>
</html>"));
      end linearPotentialDiff;
    end Kinetics;

    partial model PartialSpeciationWithSubunitsDefinition "Chemical speciation"
      extends Chemical.Processes.Internal.PartialSpeciation;

      import Chemical.Substances.Liquid;
      import Chemical.Substances.Gas;
      import Chemical.Substances.Solid;
      import Chemical.Substances.Aqueous;

      import FirstProductChoice =
             Chemical.Utilities.Types.FirstProductChoice;

      parameter FirstProductChoice firstForeSubunitFrom=FirstProductChoice.Process  "First product definition comes from?"
          annotation(HideResult=true, Dialog(group="Conditional inputs"));

      parameter Chemical.Interfaces.Definition firstForeSubunit=dropOfCommons.DefaultSubstance "First product definition"
        annotation (choicesAllMatching=true, Dialog(enable=(firstForeSubunitFrom == FirstProductChoice.Substance)));

      // for unknown reason this array can not be empty in Dymola 2025 - stupid fix: max(1,nP-1) instead of max(0,nP-1) is used as size of the array
      parameter Chemical.Interfaces.Definition nextForeSubunits[:]=fill(dropOfCommons.DefaultSubstance, max(1,nP-1)) "Definitions of next subunit_fore"
        annotation (choicesAllMatching=true, Dialog(enable=(nP > 1)));

      parameter Chemical.Interfaces.Definition process = Chemical.Interfaces.processData(1)
       "Process changes of Gibbs energy, enthalpy, volume and heat capacity (subunit_fore - reactants)"
       annotation (Dialog(enable=(firstForeSubunitFrom == FirstProductChoice.Process)));

    /*
  parameter Real K=1 "Process dissociation constant at 25°C,1bar"
   annotation (HideResult=true, Dialog(group="Process definition",enable=(firstForeSubunitFrom == FirstProductChoice.ProcessProperties)));

  parameter Modelica.Units.SI.MolarEnergy dH=0 "Process molar enthalpy change at 25°C,1bar"
   annotation (HideResult=true,Dialog(group="Process definition",enable=(firstForeSubunitFrom == FirstProductChoice.ProcessProperties)));

  parameter Modelica.Units.SI.MolarHeatCapacity dCp=0 "Process molar heat capacity change at 25°C,1bar"
   annotation (HideResult=true,Dialog(group="Process definition",enable=(firstForeSubunitFrom == FirstProductChoice.ProcessProperties)));

  parameter Modelica.Units.SI.SpecificVolume dVs=0 "Process specific volume change at 25°C,1bar [L/g]"
   annotation (HideResult=true,Dialog(group="Process definition",enable=(firstForeSubunitFrom == FirstProductChoice.ProcessProperties)));
*/
    equation

      if (nP>1) then
        for i in 1:nP-1 loop
           subunit_fore[i+1].definition = nextForeSubunits[i];
        end for;
      end if;

      if (nP>0)  then

        if (firstForeSubunitFrom == FirstProductChoice.Substance)  then
          subunit_fore[1].definition = firstForeSubunit;
        elseif (nP>1) then
          for i in 1:1 loop //this stupid loop is only for disabling check error about empty array p[:] in Dymola 2025
            subunit_fore[1].definition =
               (ones(nS)*subunit_rear.definition + process - ones(max(1,nP-1))*nextForeSubunits);
          end for;

        else
          for i in 1:1 loop //this stupid loop is only for disabling check error about empty array p[:] in Dymola 2025
            subunit_fore[1].definition = (ones(nS)*subunit_rear.definition + process);
          end for;
        end if;
      end if;

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}})),
        Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichiometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
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
<td><p>stoichiometric coefficients of i-th substance</p></td>
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
    end PartialSpeciationWithSubunitsDefinition;

    partial model PartialSpeciation "Chemical Speciation"

      import Chemical.Utilities.Types.SolutionChoice;
      import Chemical;
      import Chemical.Utilities.Types.InitializationMethods;

      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
        annotation(HideResult=true, Dialog(tab="Advanced"));

      parameter Modelica.Units.SI.MolarFlowRate n_flow_coef_reg=dropOfCommons.n_flow_coef_reg "Regulation forward/backward flow tolerance of chemical kinectics - smallest significant rate near forward is (1-coef_reg)*q_f (resp. (1-coef_reg)*q_b for backward)";
      parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
        annotation(HideResult=true, Dialog(tab="Advanced"));
      parameter InitializationMethods initN_flow =Chemical.Utilities.Types.InitializationMethods.none  "Initialization method for n_flow"
        annotation(HideResult=true, Dialog(tab= "Initialization", group="Molar flow"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_0 = 0 "Initial value for n_flow"
        annotation(HideResult=true, Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.state)));
      parameter Utilities.Units.MolarFlowAcceleration n_acceleration_0 = 0 "Initial value for der(n_flow)"
        annotation(HideResult=true, Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.derivative)));
      parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (HideResult=true, Dialog(tab="Advanced"));
      parameter Utilities.Units.Inertance L = dropOfCommons.L "Inertance of the flow"
        annotation(HideResult=true, Dialog(tab="Advanced"));
      Modelica.Units.SI.MolarFlowRate rr(stateSelect=n_flowStateSelect) "Molar flow rate of change macromolecule from rear state to fore state";
      parameter Integer nS=0 "Number of subunits in rear macromolecule state"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));
      parameter Integer nP=0 "Number of subunits in fore macromolecule state"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));
      Chemical.Interfaces.Rear subunit_rear[nS] "Subunits in rear macromolecule state" annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-100,0}), iconTransformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-100,0})));
      Chemical.Interfaces.Fore subunit_fore[nP] "Subunits in fore macromolecule state" annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={100,0}), iconTransformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={100,0})));
      Modelica.Units.SI.TemperatureSlope dT "Temperature change";
      Interfaces.SolutionWhole solution_rear(
         dT=dT,
         T=solutionState.T,
         p=solutionState.p,
         v=solutionState.v,
         n=solutionState.n,
         m=solutionState.m,
         V=solutionState.V,
         G=solutionState.G,
         Q=solutionState.Q,
         I=solutionState.I)
         "To connect solution of subunits with rear state of macromolecule"
        annotation (Placement(transformation(extent={{-110,34},{-90,54}}),   iconTransformation(extent={{-110,36},{-90,56}})));
      Interfaces.SolutionWhole solution_fore(
         dT=dT,
         T=solutionState.T,
         p=solutionState.p,
         v=solutionState.v,
         n=solutionState.n,
         m=solutionState.m,
         V=solutionState.V,
         G=solutionState.G,
         Q=solutionState.Q,
         I=solutionState.I)
          "To connect solution of subunits with fore state of macromolecule"
        annotation (Placement(transformation(extent={{90,38},{110,58}}),     iconTransformation(extent={{90,38},{110,58}})));
      Modelica.Units.SI.AmountOfSubstance nm_fore
        "Amount of the macromolecule in fore state";
      Modelica.Units.SI.AmountOfSubstance nm_rear
        "Amount of the macromolecule in rear state";
      Modelica.Units.SI.MoleFraction xm_fore
        "Mole fraction of the macromolecule in fore state";
      Modelica.Units.SI.MoleFraction xm_rear
        "Mole fraction of the macromolecule in rear state";
      Chemical.Interfaces.SolutionPart solution(
          dT=dT,
          p=solutionState.p,
          v=solutionState.v,
          n=solutionState.n,
          m=solutionState.m,
          V=solutionState.V,
          G=solutionState.G,
          Q=solutionState.Q,
          I=solutionState.I,
          i=i,
          dH=dH,
          dV=dV,
          nj=nj,
          mj=mj,
          Vj=Vj,
          Gj=Gj,
          Qj=Qj,
          Ij=Ij)
             "To connect substance with solution, where is present"
        annotation (Placement(transformation(extent={{-66,-110},{-46,-90}}), iconTransformation(extent={{-66,-110},{-46,-90}})));
         Chemical.Interfaces.SolutionState solutionState;
         Chemical.Interfaces.Definition macromolecule_rear "Specific rear state of macromolecule composed only with connected subunits";
         Chemical.Interfaces.Definition macromolecule_fore "Specific fore state of macromolecule composed only with connected subunits";
         Real i,dH,dV,nj,mj,Vj,Gj,Qj,Ij;
         Real RTlnxm_rear[nS],RTlnxm_fore[nP];
    protected
      Modelica.Units.SI.MolarEnthalpy h_fore_mix, h_rear_mix;
      Real Sx_fore,Px_rear,Kx;
      Modelica.Units.SI.ChemicalPotential uPure_subunit_rear[nS];
      Modelica.Units.SI.ChemicalPotential uPure_subunit_fore[nP];
      outer DropOfCommons dropOfCommons;
      Real Px_fore,Sx_rear;
      Real _rR[nS+nP];
      Real _uIn[nS+nP];
      Real _uOut[nS+nP];
      Real _qIn[nS+nP];
      Real _uS[nS],_uP[nP];
      Real du;
    initial equation
      if initN_flow == InitializationMethods.state then
        rr = n_flow_0;
      elseif initN_flow == InitializationMethods.derivative then
        der(rr) = n_acceleration_0;
      elseif initN_flow == InitializationMethods.steadyState then
        der(rr) = 0;
      end if;
    equation
      solutionState.T = solution.T;
      if nS>0 then
        subunit_fore.solution_forwards = fill(solutionState,nP);
        subunit_rear.solution_rearwards = fill(solutionState,nS);
      end if;
      solution_rear.dH + solution_fore.dH + dH = 0 "enthalpy change per each subunit";
      solution_rear.i + solution_fore.i + i = 0 "current change per each subunit";
      solution_rear.Qj + solution_fore.Qj + Qj = 0 "electric charge per each subunit";
      solution_rear.Ij + solution_fore.Ij + Ij = 0 "ionic strength per each subunit";
      solution_rear.nj/max(Modelica.Constants.eps,nS) + solution_fore.nj/max(Modelica.Constants.eps,nP) + nj = 0 "amount of part of particle per amount of particles";
      solution_rear.mj + solution_fore.mj + mj = 0 "mass per each subunit";
      solution_rear.Vj + solution_fore.Vj + Vj = 0 "volume per each subunit";
      solution_rear.Gj + solution_fore.Gj + Gj = 0 "Gibbs energy per each subunit";
      solution_rear.dV + solution_fore.dV + dV = 0 "volume change per each subunit";
      nm_rear*max(Modelica.Constants.eps,nS) + solution_rear.nj = 0;
      nm_fore*max(Modelica.Constants.eps,nP) + solution_fore.nj = 0;
      xm_rear = nm_rear/solutionState.n;
      xm_fore = nm_fore/solutionState.n;
      macromolecule_rear = ones(nS) * subunit_rear.definition;
      macromolecule_fore = ones(nP) * subunit_fore.definition;
      RTlnxm_rear = Modelica.Constants.R*solution.T*log(xm_rear)*ones(nS);
      RTlnxm_fore = Modelica.Constants.R*solution.T*log(xm_fore)*ones(nP);
      Sx_fore = xm_rear;
      Px_rear = xm_fore;
      Sx_rear = xm_rear;
      Px_fore = xm_fore;
      Kx = 1;
      rr*ones(nS) = subunit_rear.n_flow;
      rr*ones(nP) = -subunit_fore.n_flow;
      subunit_fore.state_forwards.h = h_fore_mix*ones(nP);
      subunit_rear.state_rearwards.h = h_rear_mix*ones(nS);
      if nS>0 then
        h_rear_mix*(subunit_rear.n_flow*ones(nS)) + subunit_fore.n_flow*subunit_fore.state_rearwards.h = 0;
      else
        h_rear_mix = 0;
      end if;
      if nP>0 then
        h_fore_mix*(subunit_fore.n_flow*ones(nP)) + subunit_rear.n_flow*subunit_rear.state_forwards.h = 0;
      else
        h_fore_mix = 0;
      end if;
      (der(subunit_rear.n_flow)*L) =  subunit_rear.r - _rR[1:nS];
      (der(subunit_fore.n_flow)*L) =  subunit_fore.r - _rR[nS+1:end];
      for i in 1:nS loop
        _uIn[i] = subunit_rear[i].state_forwards.u;
        _uOut[i] = (_uIn*_qIn - _uIn[i]*_qIn[i])/(_qIn*ones(nS+nP)-_qIn[i]);
        _qIn[i] = max(subunit_rear[i].n_flow,n_flow_reg);
        Chemical.Utilities.Internal.regStep(subunit_rear[i].n_flow,_uIn[i],_uOut[i],n_flow_reg) + _rR[i] = _uS[i];
        subunit_rear[i].state_rearwards.u = _uOut[i];
        uPure_subunit_rear[i] = Chemical.Interfaces.Properties.electroChemicalPotentialPure(subunit_rear[i].definition, subunit_rear[i].solution_forwards);
      end for;
      for i in 1:nP loop
        _uIn[i+nS] = subunit_fore[i].state_rearwards.u;
        _uOut[i+nS] = (_uIn*_qIn - _uIn[i+nS]*_qIn[i+nS])/(_qIn*ones(nS+nP)-_qIn[i+nS]);
        _qIn[i+nS] = max(subunit_fore[i].n_flow,n_flow_reg);
        Chemical.Utilities.Internal.regStep(subunit_fore[i].n_flow,_uIn[i+nS],_uOut[i+nS],n_flow_reg) + _rR[i+nS] = _uP[i];
        subunit_fore[i].state_forwards.u = _uOut[i+nS];
        uPure_subunit_fore[i] = Chemical.Interfaces.Properties.electroChemicalPotentialPure( subunit_fore[i].definition,  subunit_fore[i].solution_rearwards);
      end for;
      du = ones(nP)*_uP - ones(nS)*_uS;
        annotation(HideResult=true, Dialog(tab="Advanced"),
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}})),
        Documentation(revisions="<html>
<p><i>2013-2025 by </i>Marek Mateják </p>
</html>",     info="<html>
<h4><span style=\"color: #008000\">Notations</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>A<sub>i</sub></p></td>
<td><p>i-th substance</p></td>
</tr>
<tr>
<td><p>v<sub>i</sub></p></td>
<td><p>stoichiometric coefficients of i-th substance</p></td>
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

    end PartialSpeciation;

  end Internal;

  package Tests "Tests for top level components of undirected"
    extends Modelica.Icons.ExamplesPackage;

    model TestFlow "Test for the undirected flow resistance"
      extends Modelica.Icons.Example;

      Chemical.Boundaries.BoundaryRear boundary_rear(u0_par=100000)
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-30,0})));
      Chemical.Boundaries.BoundaryFore boundary_fore(usePotential=true, u0_par=110000) annotation (Placement(transformation(extent={{20,-10},{40,10}})));
      inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=0.01) annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      Modelica.Blocks.Sources.Step step(
        height=-80000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{60,-6},{48,6}})));
      Chemical.Boundaries.BoundaryRear boundary_rear1(usePotential=true)
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-28,-38})));
      Chemical.Boundaries.BoundaryFore boundary_fore1(usePotential=false, u0_par=100000)
        annotation (Placement(transformation(extent={{22,-48},{42,-28}})));
      Reaction reaction(
        firstProductFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Reaction reaction1(
        firstProductFrom=Chemical.Utilities.Types.FirstProductChoice.Substance,
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
    end TestFlow;

    model Diffusion

      extends Modelica.Icons.Example;

      Chemical.Boundaries.Substance s2(useFore=true, mass_start=0.6) annotation (Placement(transformation(extent={{-174,20},{-154,40}})));
      Chemical.Boundaries.Substance p2(useRear=true, mass_start=0.4) annotation (Placement(transformation(extent={{-66,20},{-46,40}})));
      Chemical.Processes.Diffusion d2(initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
                                      solutionFrom=Chemical.Utilities.Types.SolutionChoice.FirstSubstrate) annotation (Placement(transformation(extent={{-118,20},{-98,40}})));
      Chemical.Boundaries.Substance s1(useFore=true, mass_start=0.6) annotation (Placement(transformation(extent={{-174,58},{-154,78}})));
      Chemical.Boundaries.Substance p1(useRear=true, mass_start=0.4) annotation (Placement(transformation(extent={{-66,58},{-46,78}})));
      Chemical.Processes.Diffusion d1(initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
                                      solutionFrom=Chemical.Utilities.Types.SolutionChoice.Parameter) annotation (Placement(transformation(extent={{-118,58},{-98,78}})));
      Chemical.Boundaries.Substance s3(
        useFore=true,
        useSolution=true,
        mass_start=0.6) annotation (Placement(transformation(extent={{-170,-54},{-150,-34}})));
      Chemical.Boundaries.Substance p3(
        useRear=true,
        useSolution=true,
        mass_start=0.4) annotation (Placement(transformation(extent={{-62,-54},{-42,-34}})));
      Chemical.Processes.Diffusion d3(initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
                                      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort) annotation (Placement(transformation(extent={{-114,-56},{-94,-36}})));
      Solution solution annotation (Placement(transformation(extent={{-222,-122},{-14,-12}})));
      inner Modelica.Fluid.System system annotation (Placement(transformation(extent={{-210,64},{-190,84}})));
      Chemical.Boundaries.Substance s4(
        useFore=true,
        useSolution=false,
        mass_start=0.6) annotation (Placement(transformation(extent={{78,64},{98,84}})));
      Chemical.Boundaries.Substance p4(
        useRear=true,
        useSolution=false,
        mass_start=0.4) annotation (Placement(transformation(extent={{186,64},{206,84}})));
      Chemical.Processes.Diffusion d4(initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
                                      solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort) annotation (Placement(transformation(extent={{134,62},{154,82}})));
      Solution solution1 annotation (Placement(transformation(extent={{26,-4},{234,106}})));
      Chemical.Boundaries.Substance solvent(useFore=false, useSolution=true,
        mass_start=1)                                                        annotation (Placement(transformation(extent={{194,24},{214,44}})));
      Chemical.Boundaries.Substance s5(
        useFore=true,
        useSolution=true,
        mass_start=0.6) annotation (Placement(transformation(extent={{64,-54},{84,-34}})));
      Chemical.Boundaries.Substance p5(
        useRear=true,
        useSolution=true,
        mass_start=0.4) annotation (Placement(transformation(extent={{172,-54},{192,-34}})));
      Chemical.Processes.Diffusion d5(initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
                                      solutionFrom=Chemical.Utilities.Types.SolutionChoice.FirstSubstrate) annotation (Placement(transformation(extent={{120,-56},{140,-36}})));
      Solution solution2 annotation (Placement(transformation(extent={{12,-122},{220,-12}})));
    equation
      connect(s2.fore, d2.rear) annotation (Line(
          points={{-154,30},{-118,30}},
          color={158,66,200},
          thickness=0.5));
      connect(d2.fore, p2.rear) annotation (Line(
          points={{-98,30},{-66,30}},
          color={158,66,200},
          thickness=0.5));
      connect(s1.fore, d1.rear) annotation (Line(
          points={{-154,68},{-118,68}},
          color={158,66,200},
          thickness=0.5));
      connect(d1.fore, p1.rear) annotation (Line(
          points={{-98,68},{-66,68}},
          color={158,66,200},
          thickness=0.5));
      connect(s3.fore, d3.rear) annotation (Line(
          points={{-150,-44},{-122,-44},{-122,-46},{-114,-46}},
          color={158,66,200},
          thickness=0.5));
      connect(d3.fore, p3.rear) annotation (Line(
          points={{-94,-46},{-70,-46},{-70,-44},{-62,-44}},
          color={158,66,200},
          thickness=0.5));
      connect(solution.solution, d3.solution) annotation (Line(points={{-55.6,-120.9},{-55.6,-126},{-110,-126},{-110,-56}}, color={127,127,0}));
      connect(s3.solution, solution.solution) annotation (Line(points={{-166,-54},{-166,-120.9},{-55.6,-120.9}}, color={127,127,0}));
      connect(p3.solution, solution.solution) annotation (Line(points={{-58,-54},{-58,-126},{-55.6,-126},{-55.6,-120.9}}, color={127,127,0}));
      connect(s4.fore, d4.rear) annotation (Line(
          points={{98,74},{126,74},{126,72},{134,72}},
          color={158,66,200},
          thickness=0.5));
      connect(d4.fore, p4.rear) annotation (Line(
          points={{154,72},{178,72},{178,74},{186,74}},
          color={158,66,200},
          thickness=0.5));
      connect(solution1.solution, d4.solution) annotation (Line(points={{192.4,-2.9},{192.4,-8},{138,-8},{138,62}}, color={127,127,0}));
      connect(solution1.solution, solvent.solution) annotation (Line(points={{192.4,-2.9},{192.4,-10},{198,-10},{198,24}},     color={127,127,0}));
      connect(s5.fore, d5.rear) annotation (Line(
          points={{84,-44},{112,-44},{112,-46},{120,-46}},
          color={158,66,200},
          thickness=0.5));
      connect(d5.fore, p5.rear) annotation (Line(
          points={{140,-46},{164,-46},{164,-44},{172,-44}},
          color={158,66,200},
          thickness=0.5));
      connect(s5.solution, solution2.solution) annotation (Line(points={{68,-54},{68,-120.9},{178.4,-120.9}}, color={127,127,0}));
      connect(p5.solution, solution2.solution) annotation (Line(points={{176,-54},{176,-126},{178.4,-126},{178.4,-120.9}}, color={127,127,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-240,-140},{240,140}})), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-240,
                -140},{240,140}})),
        experiment(StopTime=1, __Dymola_Algorithm="Dassl"));
    end Diffusion;

    model SimpleFlow
      extends Modelica.Icons.Example;
      Chemical.Boundaries.Substance substance(
        useFore=true,
        preferMass=false,
        amountOfSubstance_start=2)   annotation (Placement(transformation(extent={{-70,14},{-50,34}})));
      Chemical.Boundaries.Substance substance1(
        useRear=true,
        preferMass=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{48,12},{68,32}})));
      inner DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{-76,66},{-56,86}})));
      Pump pump(solutionFrom=Chemical.Utilities.Types.SolutionChoice.FirstSubstrate, SubstanceFlow=1)
        annotation (Placement(transformation(extent={{-14,14},{6,34}})));
    equation
      connect(substance.fore, pump.rear) annotation (Line(
          points={{-50,24},{-14,24}},
          color={158,66,200},
          thickness=0.5));
      connect(pump.fore, substance1.rear) annotation (Line(
          points={{6,24},{40,24},{40,22},{48,22}},
          color={158,66,200},
          thickness=0.5));
       annotation (experiment(StopTime=1, __Dymola_Algorithm="Dassl"));
    end SimpleFlow;

    model TestMembrane

      extends Modelica.Icons.Example;

      Chemical.Boundaries.Substance s1(useFore=true, mass_start=0.6) annotation (Placement(transformation(extent={{-174,58},{-154,78}})));
      Chemical.Boundaries.Substance p1(useRear=true, mass_start=0.4) annotation (Placement(transformation(extent={{-66,58},{-46,78}})));
      Membrane d1(solutionFrom=Chemical.Utilities.Types.SolutionChoice.Parameter) annotation (Placement(transformation(extent={{-118,58},{-98,78}})));
    equation
      connect(d1.fore, p1.rear) annotation (Line(
          points={{-98,68},{-66,68}},
          color={158,66,200},
          thickness=0.5));
      connect(s1.fore, d1.rear) annotation (Line(
          points={{-154,68},{-118,68}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-240,-140},{240,140}})), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-240,
                -140},{240,140}})),
        experiment(StopTime=1, __Dymola_Algorithm="Dassl"));
    end TestMembrane;
    annotation (Documentation(info="<html>
<u>Tests for top level components of the undirected chemical simulation package.</u>
</html>"));
  end Tests;

  annotation (Documentation(info="<html>
<u>This package contains models implementing undirected versions of the processes. Here, the thermodynamic state of one or more fluid streams is changed by exchanging heat or work with the streams, or by delaying the state.</u>
</html>", revisions="<html>
<u><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</u>
</html>"), Icon(graphics={
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
