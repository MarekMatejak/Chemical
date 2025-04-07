within Chemical;
package Processes "Undirected process package"
  model Reaction "Chemical Reaction"
    extends Chemical.Processes.Internal.PartialReactionWithSubstanceData;
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = if nS==0 then SolutionChoice.fromParameter else SolutionChoice.fromSubstrate);
    extends Chemical.Interfaces.ConditionalKinetics
                                          (k_forward=1);


    import Chemical.Utilities.Types.SolutionChoice;

    Real rr_fore_exact,rr_rear_exact,  kb;

    replaceable function uLoss =
        Internal.Kinetics.pleaseSelectPotentialLoss
      constrainedby
        Internal.Kinetics.partialPotentialLoss "Electro-chemical potential loss function"
      annotation(choicesAllMatching=true, Documentation(info="<html>
    <p>Electro-chemical potential loss function used in the diffusion.</p>
    </html>"));

    Real du_fore_,du_rear_;
  equation

    connect(substrates[1].solution,inputSubstrateSolution);


    du_fore_ = -uLoss(rr,Sx_fore,solutionState);
  //  if nS>0 then
      du_rear_ = -uLoss(-rr,Px_rear,solutionState);
  //  end if;

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
</html>",   info="<html>
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
</html>",   info="<html>
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

  model FastReaction "Chemical Reaction"
    extends Chemical.Processes.Internal.PartialReactionWithSubstanceData;

    parameter Real kC = 1;


  equation

    rr = kC * du_fore;
    if nS>0 then
      rr = kC * du_rear;
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
  end FastReaction;

  model FastReactionWithSolutions "Chemical Reaction"
    extends Internal.PartialReaction;

    parameter stateOfMatter.SubstanceDataParameters productsSubstanceData[nP]
     annotation (choicesAllMatching = true);

    parameter Real kC = 1;

    Interfaces.SolutionPort productSolution[nP] annotation (Placement(transformation(extent={{90,-50},{110,-30}}), iconTransformation(extent={{90,-50},{110,-30}})));
  equation

    rr = kC * du_fore;
    if nS>0 then
      rr = kC * du_rear;
    end if;

    products.definition = productsSubstanceData;

    for i in 1:nP loop
      products[i].solution.T = productSolution[i].T;
      products[i].solution.p = productSolution[i].p;
      products[i].solution.v = productSolution[i].v;
      products[i].solution.n = productSolution[i].n;
      products[i].solution.m = productSolution[i].m;
      products[i].solution.V = productSolution[i].V;
      products[i].solution.G = productSolution[i].G;
      products[i].solution.Q = productSolution[i].Q;
      products[i].solution.I = productSolution[i].I;

      productSolution[i].dH = 0;
      productSolution[i].i = 0;
      productSolution[i].Qj = 0;
      productSolution[i].Ij = 0;
      productSolution[i].nj = 0;
      productSolution[i].mj = 0;
      productSolution[i].Vj = 0;
      productSolution[i].Gj = 0;
      productSolution[i].dV = 0;
    end for;

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
</html>",   info="<html>
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
  end FastReactionWithSolutions;

  model Diffusion "Solute diffusion"
    extends Icons.Diffusion;
    extends Chemical.Interfaces.SISO(redeclare package stateOfMatterRear = stateOfMatter, redeclare package stateOfMatterFore = stateOfMatter);
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = SolutionChoice.fromParameter);

    import Chemical.Utilities.Types.SolutionChoice;

    replaceable package stateOfMatter = Interfaces.Incompressible constrainedby
      Interfaces.StateOfMatter "Substance model of inlet"
      annotation (choices(
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.Incompressible  "Incompressible"),
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

    replaceable function uLoss =
        Internal.Kinetics.pleaseSelectPotentialLoss
      constrainedby
        Internal.Kinetics.partialPotentialLoss "Electro-chemical potential loss function"
      annotation(choicesAllMatching=true, Documentation(info="<html>
    <p>Electro-chemical potential loss function used in the diffusion.</p>
    </html>"));

  equation
    fore.definition = rear.definition;
    fore.solution = solutionState;

    connect(rear.solution,inputSubstrateSolution);


    du_rear = -uLoss(n_flow,x_fore,rear.solution);
    du_fore = -uLoss(-n_flow,x_fore,fore.solution);

     annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
  end Diffusion;

  model GasSolubility "Henry's law of gas solubility into liquid."

    extends Icons.GasSolubility;
    extends Interfaces.PartialGasToLiquid;
    extends Interfaces.ConditionalKinetics(k_forward=1);

  equation

    n_flow = kf * x_in * ( 1 -  exp(-duRT));

    annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
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

  model GasVolatility "Henry's law of gas volatility from liquid."

    extends Icons.GasSolubility;
    extends Interfaces.PartialLiquidToGas;
    extends Interfaces.ConditionalKinetics(k_forward=1);

  equation

    n_flow = kf * x_in * ( 1 -  exp(-duRT));

    annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
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
  end GasVolatility;

  model FastGasSolubility "Henry's law of gas solubility in liquid."

    extends Icons.GasSolubility;
    extends Interfaces.PartialGasToLiquid;

    parameter Real kC = 1;

  equation

    n_flow = kC * du;

    annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
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
  end FastGasSolubility;

  model FastGasVolatility "Henry's law of gas solubility in liquid."

    extends Icons.GasSolubility;
    extends Interfaces.PartialLiquidToGas;

    parameter Real kC = 1;

  equation

    n_flow = kC * du;

    annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
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
  end FastGasVolatility;

  model Membrane "Passive transport of the substance through semipermeable membrane"
    extends Icons.Membrane;
    extends Interfaces.PartialChangeSolution;
    extends Interfaces.ConditionalKinetics(k_forward=1);

  equation
    //the main equation

    n_flow = kf * x_in * ( 1 -  exp(-duRT));

    annotation ( Documentation(info="<html>
<p><u><b><font style=\"color: #008000; \">Filtration throught semipermeable membrane.</font></b></u></p>
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
            lineColor={128,0,255},
          origin={69,2},
          rotation=90)}));
  end Membrane;

  model Pump "Prescribed sunstance molar flow"
    extends Interfaces.PartialChangeSolution;
    extends Interfaces.ConditionalSubstanceFlow;

  equation
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
            lineColor={128,0,255},
            origin={0,-72},
            rotation=360,
            textString="%name")}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end Pump;

  model Stream "Flow of whole solution"
    extends Interfaces.PartialChangeSolution;
    extends Onedirectional.Boundaries.Internal.ConditionalSolutionFlow;

  equation
    n_flow = c_in * volumeFlow;

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

  model SpeciationIn "Quaternary macromolecule form defined by all its subunits"
    extends Icons.Speciation;

    replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby
      Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
       annotation (choicesAllMatching = true);

    /*parameter stateOfMatter.SubstanceDataParameters substanceData
  "Definition of the substance"
    annotation (choicesAllMatching = true);
*/
    parameter stateOfMatter.SubstanceDataParameters subunitsSubstanceData[NumberOfSubunits]
    "Definition of the subunits"
     annotation (choicesAllMatching = true);

    parameter Integer NumberOfSubunits=1
    "Number of independent subunits occurring in macromolecule";

    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Utilities.Units.Inertance L = dropOfCommons.L "Inertance of the flow"
      annotation(Dialog(tab="Advanced"));

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
    Interfaces.Inlet         inlet annotation (Placement(transformation(extent={{110,-110},{90,-90}}), iconTransformation(extent={{110,-110},{90,-90}})));
    Interfaces.Outlet subunits[NumberOfSubunits] "Subunits of macromolecule" annotation (Placement(transformation(extent={{-56,-14},{-36,66}}),
          iconTransformation(
          extent={{-10,-40},{10,40}},
          rotation=90,
          origin={-30,102})));

    parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));

  protected
    outer DropOfCommons dropOfCommons;
    Chemical.Utilities.Units.URT uRT_out;
    Modelica.Units.SI.ChemicalPotential u_out;
  equation

    inlet.r + inlet.state.u = u_out;

    if NumberOfSubunits>0 then
      (ones(NumberOfSubunits) * subunits.r) = (inlet.r)  -  der(inlet.n_flow)*L;
      for i in 2:NumberOfSubunits loop
        //first subunit is based on inertial potential,
        //other subunits are provided as source
        der(subunits[i].state.u).*TC = subunits[i].r;
      end for;
    end if;

    //amount of macromolecule (all forms in conformation)
    nm*NumberOfSubunits + subunitSolution.nj = 0;

    //change of macromolecule = change of its subunits
    subunits.n_flow = -inlet.n_flow * ones(NumberOfSubunits);

    //mole fraction of all forms in conformation
    xm = nm/solution.n;

    //electrochemical potential of the specific form divided by (Modelica.Constants.R*solution.T)
    uRT_out = log(xm) +
          sum(subunits.state.u/(Modelica.Constants.R*inlet.solution.T) - log(xm)
           * ones(NumberOfSubunits));

    u_out = uRT_out*(Modelica.Constants.R*inlet.solution.T);

    subunits.state.h = (inlet.state.h/NumberOfSubunits)*ones(NumberOfSubunits);

    //inlet.state.u = sum(subunits.state.u - (log(xm)*(Modelica.Constants.R*inlet.solution.T)) * ones(NumberOfSubunits));

    for i in 1:NumberOfSubunits loop
      subunits[i].solution.T = solution.T;
      subunits[i].solution.v = solution.v;
      subunits[i].solution.p = solution.p;
      subunits[i].solution.n = solution.n;
      subunits[i].solution.m = solution.m;
      subunits[i].solution.V = solution.V;
      subunits[i].solution.G = solution.G;
      subunits[i].solution.Q = solution.Q;
      subunits[i].solution.I = solution.I;
    end for;

    //inlet.definition = substanceData;
    subunits.definition = subunitsSubstanceData;

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
</html>",   info="<html>
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
  end SpeciationIn;

  model SpeciationOut "Quaternary macromolecule form defined by all its subunits"
    extends Icons.Speciation;

    replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby
      Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
       annotation (choicesAllMatching = true);

    parameter stateOfMatter.SubstanceDataParameters substanceData
    "Definition of the substance"
      annotation (choicesAllMatching = true);

    parameter Integer NumberOfSubunits=1
    "Number of independent subunits occurring in macromolecule";

    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Utilities.Units.Inertance L = dropOfCommons.L "Inertance of the flow"
      annotation(Dialog(tab="Advanced"));

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
    Interfaces.Outlet outletSubstance annotation (Placement(transformation(extent={{90,-110},{110,-90}}), iconTransformation(extent={{90,-110},{110,-90}})));
    Interfaces.Inlet subunits[NumberOfSubunits] "Subunits of macromolecule" annotation (Placement(transformation(extent={{-56,-14},{-36,66}}),
          iconTransformation(
          extent={{10,-40},{-10,40}},
          rotation=90,
          origin={-30,102})));

    parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));

  //protected

   // stateOfMatter.SubstanceData outletSubstanceDefinition;
    outer DropOfCommons dropOfCommons;
  //  Modelica.Units.SI.MolarEnthalpy h_out;
    Chemical.Utilities.Units.URT uRT_out;
    Modelica.Units.SI.ChemicalPotential u_out; //, u_pure;
  //  Modelica.Units.SI.MolarEntropy s_pure;

  equation

    outletSubstance.state.u = u_out;

    if NumberOfSubunits>0 then
      (ones(NumberOfSubunits) * subunits.r) =(outletSubstance.r) - der(outletSubstance.n_flow)*L;
    end if;

    //amount of macromolecule (all forms in conformation)
    nm*NumberOfSubunits + subunitSolution.nj = 0;

    //change of macromolecule = change of its subunits
    subunits.n_flow =-outletSubstance.n_flow*ones(NumberOfSubunits);

    //mole fraction of all forms in conformation
    xm = nm/solution.n;

    //electrochemical potential of the specific form divided by (Modelica.Constants.R*solution.T)
    uRT_out = log(xm) +
          sum(subunits.state.u/(Modelica.Constants.R*outletSubstance.solution.T) - log(xm)
           * ones(NumberOfSubunits));

    u_out = uRT_out*(Modelica.Constants.R*outletSubstance.solution.T);

    subunits.state.h = (outletSubstance.state.h/NumberOfSubunits)*ones(NumberOfSubunits);

   // h_out =outletSubstance.state.h;
   //  (subunits.state.h*ones(NumberOfSubunits)) =(outletSubstance.state.h);

    outletSubstance.solution.T = solution.T;
    outletSubstance.solution.v = solution.v;
    outletSubstance.solution.p = solution.p;
    outletSubstance.solution.n = solution.n;
    outletSubstance.solution.m = solution.m;
    outletSubstance.solution.V = solution.V;
    outletSubstance.solution.G = solution.G;
    outletSubstance.solution.Q = solution.Q;
    outletSubstance.solution.I = solution.I;

    outletSubstance.definition = substanceData;

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
</html>",   info="<html>
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
  end SpeciationOut;

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
    partial model PartialReactionWithSubstanceData "Chemical Reaction"
      extends Chemical.Processes.Internal.PartialReaction;

      parameter stateOfMatter.SubstanceDataParameters productsSubstanceData[nP] "Array of products definitions"
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

      Chemical.Interfaces.Rear substrates[nS](redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-100,0}), iconTransformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-100,0})));

      Chemical.Interfaces.Fore products[nP](redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(
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
</html>",     info="<html>
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

    partial model PartialSwitchSolution "Substance between different chemical solutions"
      extends Chemical.Interfaces.SISO;

      Interfaces.SolutionPort solution annotation (Placement(transformation(extent={{92,-52},{112,-32}}),   iconTransformation(extent={{92,-52},{112,-32}})));

    equation

      fore.definition = rear.definition;

      fore.solution.T = solution.T;
      fore.solution.p = solution.p;
      fore.solution.v = solution.v;
      fore.solution.n = solution.n;
      fore.solution.m = solution.m;
      fore.solution.V = solution.V;
      fore.solution.G = solution.G;
      fore.solution.Q = solution.Q;
      fore.solution.I = solution.I;

      solution.dH = 0;
      solution.i = 0;
      solution.Qj = 0;
      solution.Ij = 0;
      solution.nj = 0;
      solution.mj = 0;
      solution.Vj = 0;
      solution.Gj = 0;
      solution.dV = 0;

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}})),
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
    end PartialSwitchSolution;

    package Kinetics
      function pleaseSelectPotentialLoss "Please select chemical potential loss function"
        extends partialPotentialLoss;

      algorithm
        assert(false, "Please select potential loss function");

        potentialLoss :=0;

        annotation (Documentation(info="<html>
<p>
Potential loss function without actual equations with an always failing assert to
output a meaningful error, when the user forgot to select a function. This should
be used as a default.
</p>
</html>"));
      end pleaseSelectPotentialLoss;

      partial function partialPotentialLoss "Interface for potential loss functions"
        extends Modelica.Icons.Function;

        input Modelica.Units.SI.MolarFlowRate n_flow "Molar flow rate";
        input Modelica.Units.SI.MoleFraction x "Substrates mole fraction product";
        input Chemical.Interfaces.SolutionState solutionState "Solution state";

        output Modelica.Units.SI.ChemicalPotential potentialLoss "Gibbs energy lost in chemical process";

        annotation(Inline=true, smoothOrder=100,
          Documentation(info="<html>
    <p>Interface definition for a potential loss in a chemical process. Inputs are information about flow condition and the chemical solution state, output is the electro-chemical potential drop.</p>
</html>"));
      end partialPotentialLoss;

      function classicPotentialLoss "Classical potential loss function"
        extends partialPotentialLoss;

        input Real kf(unit="Pa.s/kg") = 1 "Forward rate coeeficient"
          annotation(Dialog(enable=true));

      algorithm
        potentialLoss := Modelica.Constants.R*solutionState.T*log(1-n_flow/(kf*x));

        annotation (Documentation(info="<html>
<p>
This Gibbs energy loss du to reach chemical process molar flow q: 
</p>
<blockquote><pre>
 = kf * xA - kb * xB
</pre></blockquote>
<p>
where kf is the forwar rate coefficient and kb is backward rate coefficient of chemical process;
xA is mole fraction os substrate and xB is mole fraction of product.
And K = kf/kb = xB/xA is a dissociation coefficient of the chemical process.
</p>
</html>"));
      end classicPotentialLoss;

      function fastPotentialLoss "Fast potential loss function"
        extends partialPotentialLoss;

        input Real kC(unit="mol2/(J.s)") = 1 "Linear factor"
          annotation(Dialog(enable=true));


      algorithm
        potentialLoss := n_flow/kC;

        annotation (Documentation(info="<html>
<p>
This Gibbs energy loss is linear in the molarflow with the linear factor kC: 
</p>
<blockquote><pre>
du := n_flow/kC;
</pre></blockquote>
</html>"));
      end fastPotentialLoss;
    end Kinetics;
  end Internal;

  package Tests "Tests for top level components of undirected"
    extends Modelica.Icons.ExamplesPackage;

    model TestFlowResistance "Test for the undirected flow resistance"
      extends Modelica.Icons.Example;

      Chemical.Boundaries.BoundaryRear boundary_rear(substanceData=Chemical.Substances.Water_liquid(), u0_par=100000)
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-30,0})));
      Chemical.Boundaries.BoundaryFore boundary_fore(potentialFromInput=true, u0_par=110000) annotation (Placement(transformation(extent={{20,-10},{40,10}})));
      inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=0.01) annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      Modelica.Blocks.Sources.Step step(
        height=-80000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{60,-6},{48,6}})));
      Chemical.Boundaries.BoundaryRear boundary_rear1(substanceData=Chemical.Substances.Water_liquid(), potentialFromInput=true)
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-28,-38})));
      Chemical.Boundaries.BoundaryFore boundary_fore1(potentialFromInput=false, u0_par=100000)
        annotation (Placement(transformation(extent={{22,-48},{42,-28}})));
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
      import Chemical;
      import Chemical;
       extends Modelica.Icons.Example;

      constant Real K = 2 "Dissociation constant of the reaction";

      constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Boundaries.Substance A(
        useRear=false,
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-52,-8},{-32,12}})));

      Chemical.Processes.FastReaction
               fastReaction(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(K))},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
      Chemical.Boundaries.Substance B(
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
      connect(A.fore, fastReaction.substrates[1]) annotation (Line(
          points={{-32,2},{-10,2}},
          color={158,66,200},
          thickness=0.5));
      connect(fastReaction.products[1], B.rear) annotation (Line(
          points={{10,2},{42,2}},
          color={158,66,200},
          thickness=0.5));
      annotation (Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A and substance B, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that mole fraction (A.x and B.x) are always summed to 1 for the solution.</p>
</html>"),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
    end SimpleReaction;

    model SimpleReaction2 "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
      import Chemical;
      import Chemical;
      import Chemical;
       extends Modelica.Icons.Example;

      constant Real Kb(unit="kg/mol") = 2
        "Molarity based dissociation constant of the reaction with one more reactant";

      constant Real Kx(unit="1") = Kb*55.508
        "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

      constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Boundaries.Substance A(
        useFore=true,
        useSolution=true,
        substanceData(MolarWeight=1),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,4},{-14,24}})));
      Reaction reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(MolarWeight=2, DfG=-R*T_25degC*log(Kx))},
        nS=2,
        nP=1) annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Chemical.Boundaries.Substance B(
        useFore=true,
        useSolution=true,
        substanceData(MolarWeight=1),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Chemical.Boundaries.Substance C(
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
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
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
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
        experiment(StopTime=100, __Dymola_Algorithm="Dassl"));
    end SimpleReaction22;

    model ExothermicReaction "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"
      import Chemical;
      import Chemical;
      import Chemical;
      import Chemical;

       extends Modelica.Icons.Example;

      parameter Modelica.Units.SI.MolarEnergy ReactionEnthalpy=-55000;

      Chemical.Solution thermal_isolated_solution(useMechanicPorts=true, ConstantTemperature=false)
        annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));
      Chemical.Boundaries.Substance A(
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
      Reaction reaction2_2(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfH=ReactionEnthalpy)},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
      Chemical.Boundaries.Substance B(
        useRear=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{20,-60},{40,-40}})));

      Chemical.Solution solution_at_constant_temperature(useMechanicPorts=true, useThermalPort=true)
        annotation (Placement(transformation(extent={{-100,0},{98,94}})));
      Chemical.Boundaries.Substance A1(
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      Reaction reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfH=ReactionEnthalpy)},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-8,40},{12,60}})));
      Chemical.Boundaries.Substance B1(
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
</html>",   info="<html>
<p>Demonstration of exotermic reaction with perfect cooling (i.e. connected fixed temperature to the HeatPort) and thermally insulated (HetPort unconnected). See solution_(...).T</p>
</html>"),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"),
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

      Chemical.Processes.ConductionElement conductionElement(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=-1,
        V(displayUnit="l") = 0.001,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.T,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,60},{10,80}})));

      Chemical.Boundaries.BoundaryFore boundary_fore(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,60},{40,80}})));
      Chemical.Boundaries.BoundaryRear boundary_rear(
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
      Chemical.Processes.ConductionElement conductionElement1(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        V(displayUnit="l") = 0.001,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.h,
        T_0=263.15,
        h_0=1000,
        neglectChemicalPotentialChanges=false) annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      Chemical.Boundaries.BoundaryFore boundary_fore1(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,30},{40,50}})));
      Chemical.Boundaries.BoundaryRear boundary_rear1(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,30},{-20,50}})));
      Modelica.Blocks.Sources.Ramp ramp(
        height=2e2,
        duration=0.001,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,30},{-50,50}})));
      Chemical.Processes.ConductionElement conductionElement2(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=-1,
        V(displayUnit="l") = 0.001,
        enforce_global_energy_conservation=true,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.fore,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,0},{10,20}})));

      Chemical.Boundaries.BoundaryFore boundary_fore2(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,0},{40,20}})));
      Chemical.Boundaries.BoundaryRear boundary_rear2(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
      Modelica.Blocks.Sources.Step step1(
        height=2e2,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
      Chemical.Processes.ConductionElement conductionElement3(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=-1,
        V(displayUnit="l") = 0.001,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.rear,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));

      Chemical.Boundaries.BoundaryFore boundary_fore3(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,-30},{40,-10}})));
      Chemical.Boundaries.BoundaryRear boundary_rear3(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
      Modelica.Blocks.Sources.Step step2(
        height=2e2,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,-30},{-50,-10}})));
      Chemical.Processes.ConductionElement conductionElement4(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=0,
        V(displayUnit="l") = 0.001,
        enforce_global_energy_conservation=true,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.port,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

      Chemical.Boundaries.BoundaryFore boundary_fore4(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,-60},{40,-40}})));
      Chemical.Boundaries.BoundaryRear boundary_rear4(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
      Modelica.Blocks.Sources.Step step3(
        height=2e2,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,-60},{-50,-40}})));
      Chemical.Processes.ConductionElement conductionElement5(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=-1,
        V(displayUnit="l") = 0.001,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.port,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,-90},{10,-70}})));

      Chemical.Boundaries.BoundaryFore boundary_fore5(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,-90},{40,-70}})));
      Chemical.Boundaries.BoundaryRear boundary_rear5(
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

      Chemical.Processes.TransportDelay transportDelay1(
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
      Chemical.Boundaries.BoundaryFore boundary_fore(redeclare package stateOfMatter = stateOfMatter, u0_par=100000)
        annotation (Placement(transformation(extent={{70,-50},{90,-30}})));
      Chemical.Boundaries.BoundaryRear boundary_rear(
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
      import Chemical;
      import Chemical;
      import Chemical;
      import Chemical;
      import Chemical;
      extends Modelica.Icons.Example;

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
      Chemical.Boundaries.Substance P(
        useRear=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=1e-8) annotation (Placement(transformation(extent={{72,-12},{92,8}})));

      Chemical.Boundaries.Substance S(
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=100) annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));

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

      Chemical.Boundaries.Substance ES(
        useRear=true,
        useFore=true,
        useSolution=true,
        initAmount=Chemical.Utilities.Types.InitializationMethods.state,
        use_mass_start=false,
        amountOfSubstance_start=tE/2) annotation (Placement(transformation(extent={{-8,-10},{12,10}})));
      Chemical.Boundaries.Substance E(
        useRear=true,
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=tE/2) annotation (Placement(transformation(extent={{12,36},{-8,56}})));
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

      Chemical.Boundaries.Substance liquidWater(
        useSolution=true,
        substanceData=Chemical.Substances.Water_liquid(),
        use_mass_start=true,
        mass_start=1) annotation (Placement(transformation(extent={{42,-80},{62,-60}})));
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
        experiment(StopTime=100000, __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);
    end EnzymeKinetics;

    model Diffusion
      Boundaries.Substance s2(useFore=true, mass_start=0.6) annotation (Placement(transformation(extent={{-174,20},{-154,40}})));
      Boundaries.Substance p2(useRear=true, mass_start=0.4) annotation (Placement(transformation(extent={{-66,20},{-46,40}})));
      Chemical.Processes.Diffusion d2(solutionFrom=Chemical.Utilities.Types.SolutionChoice.fromSubstrate, redeclare function uLoss =
            Chemical.Processes.Internal.Kinetics.classicPotentialLoss) annotation (Placement(transformation(extent={{-118,20},{-98,40}})));
      Boundaries.Substance s1(useFore=true, mass_start=0.6) annotation (Placement(transformation(extent={{-174,58},{-154,78}})));
      Boundaries.Substance p1(useRear=true, mass_start=0.4) annotation (Placement(transformation(extent={{-66,58},{-46,78}})));
      Chemical.Processes.Diffusion d1(solutionFrom=Chemical.Utilities.Types.SolutionChoice.fromParameter, redeclare function uLoss =
            Internal.Kinetics.classicPotentialLoss) annotation (Placement(transformation(extent={{-118,58},{-98,78}})));
      Boundaries.Substance s3(
        useFore=true,
        useSolution=true,
        mass_start=0.6) annotation (Placement(transformation(extent={{-170,-54},{-150,-34}})));
      Boundaries.Substance p3(
        useRear=true,
        useSolution=true,
        mass_start=0.4) annotation (Placement(transformation(extent={{-62,-54},{-42,-34}})));
      Chemical.Processes.Diffusion d3(solutionFrom=Chemical.Utilities.Types.SolutionChoice.fromSolutionPort, redeclare function uLoss =
            Internal.Kinetics.classicPotentialLoss) annotation (Placement(transformation(extent={{-114,-56},{-94,-36}})));
      Solution solution annotation (Placement(transformation(extent={{-222,-122},{-14,-12}})));
      inner Modelica.Fluid.System system annotation (Placement(transformation(extent={{-210,64},{-190,84}})));
      Boundaries.Substance s4(
        useFore=true,
        useSolution=false,
        mass_start=0.6) annotation (Placement(transformation(extent={{78,64},{98,84}})));
      Boundaries.Substance p4(
        useRear=true,
        useSolution=false,
        mass_start=0.4) annotation (Placement(transformation(extent={{186,64},{206,84}})));
      Chemical.Processes.Diffusion d4(solutionFrom=Chemical.Utilities.Types.SolutionChoice.fromSolutionPort, redeclare function uLoss =
            Internal.Kinetics.classicPotentialLoss) annotation (Placement(transformation(extent={{134,62},{154,82}})));
      Solution solution1 annotation (Placement(transformation(extent={{26,-4},{234,106}})));
      Boundaries.Substance
                solvent(useFore=false, useSolution=true) annotation (Placement(transformation(extent={{194,24},{214,44}})));
      Boundaries.Substance s5(
        useFore=true,
        useSolution=true,
        mass_start=0.6) annotation (Placement(transformation(extent={{64,-54},{84,-34}})));
      Boundaries.Substance p5(
        useRear=true,
        useSolution=true,
        mass_start=0.4) annotation (Placement(transformation(extent={{172,-54},{192,-34}})));
      Chemical.Processes.Diffusion d5(solutionFrom=Chemical.Utilities.Types.SolutionChoice.fromSubstrate, redeclare function uLoss =
            Internal.Kinetics.classicPotentialLoss) annotation (Placement(transformation(extent={{120,-56},{140,-36}})));
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
                -140},{240,140}})));
    end Diffusion;

    model Diffusion2
      Boundaries.Substance substrateInSolution1(
        useFore=true,
        useSolution=false,
        mass_start=0.5) annotation (Placement(transformation(extent={{-22,42},{-2,62}})));
      Boundaries.Substance productInSolution1(
        useRear=true,
        useSolution=false,
        mass_start=0.5) annotation (Placement(transformation(extent={{86,40},{106,60}})));
      Chemical.Processes.Diffusion diffusion3(solutionFrom=Chemical.Utilities.Types.SolutionChoice.fromSolutionPort, redeclare function uLoss =
            Chemical.Processes.Internal.Kinetics.fastPotentialLoss) annotation (Placement(transformation(extent={{34,38},{54,58}})));
      Solution solution1 annotation (Placement(transformation(extent={{-74,-28},{134,82}})));
      Boundaries.Substance
                solvent(useFore=false, useSolution=true) annotation (Placement(transformation(extent={{94,0},{114,20}})));
    equation
      connect(substrateInSolution1.fore, diffusion3.rear)
        annotation (Line(
          points={{-2,52},{26,52},{26,48},{34,48}},
          color={158,66,200},
          thickness=0.5));
      connect(diffusion3.fore, productInSolution1.rear)
        annotation (Line(
          points={{54,48},{78,48},{78,50},{86,50}},
          color={158,66,200},
          thickness=0.5));
      connect(solution1.solution, diffusion3.solution) annotation (Line(points={{92.4,-26.9},{92.4,-32},{38,-32},{38,38}}, color={127,127,0}));
      connect(solution1.solution, solvent.solution) annotation (Line(points={{92.4,-26.9},{92.4,-34},{98,-34},{98,0}}, color={127,127,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{760,100}})), Diagram(coordinateSystem(preserveAspectRatio=false, extent
              ={{-100,-100},{760,100}})));
    end Diffusion2;
    annotation (Documentation(info="<html>
<u>Tests for top level components of the undirected chemical simulation package.</u>
</html>"));
  end Tests;

  model TransportDelay "Delay chemical state depending on fluid speed"
    extends Chemical.Interfaces.SISO                 (final cliu_u_out=
          false);

    parameter Modelica.Units.SI.Length l "Length of Delay Pipe";
    parameter Modelica.Units.SI.Radius r "Radius of Delay Pipe";
    parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimal Density" annotation (Dialog(tab="Advanced"));
    parameter Real v_min(
      min=0,
      unit="1/s")=0.01                             "Minimum nondimensional speed"
      annotation(Dialog(tab="Advanced"));
    parameter Real v_max(
      min=0,
      unit="1/s")=50                             "Maximum nondimensional speed"
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
