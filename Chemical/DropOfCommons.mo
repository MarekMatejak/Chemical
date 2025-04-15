within Chemical;
model DropOfCommons "Model for global parameters"

  parameter Chemical.Utilities.Units.Inertance L=1e-3 "Inertance of the molar flow through electro-chemical process" annotation (Dialog(tab="Advanced"));

  parameter Modelica.Units.SI.Time TC=1e-3 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));

  parameter Modelica.Units.SI.MolarFlowRate n_flow_reg = 1e-5 "Regularization threshold of molar flow rate"
    annotation(Dialog(group="Regularization"));

  parameter AssertionLevel assertionLevel = AssertionLevel.error "Global assertion level";

  parameter Chemical.Interfaces.Definition DefaultSubstance=Chemical.Substances.Liquid.Unknown "Default substance"
    annotation (choicesAllMatching=true);

  parameter Modelica.Units.SI.Mass DefaultMass=1e-3 "Default initial mass of the substance";

  parameter Modelica.Units.SI.AmountOfSubstance DefaultAmount=1e-3 "Default initial amount of substance base molecules";

  annotation (

              defaultComponentName="dropOfCommons",
    defaultComponentPrefixes="inner",
    missingInnerMessage="
Your model is using an outer \"dropOfCommons\" component but
an inner \"dropOfCommons\" component is not defined.
Use Chemical.DropOfCommons in your model
to specify system properties.",Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Ellipse(
          extent={{-80,-60},{80,-100}},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None),
        Polygon(
          points={{0,100},{16,36},{80,-32},{0,-100},{-82,-30},{-18,36},{0,100}},
          lineColor={194,138,221},
          smooth=Smooth.Bezier,
          fillColor={158,66,200},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{6,42},{20,16},{44,-14},{22,-38},{6,42}},
          smooth=Smooth.Bezier,
          fillColor={194,138,221},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None),
        Polygon(
          points={{-6,-76},{-40,-62},{-56,-30},{-30,-44},{-6,-76}},
          pattern=LinePattern.None,
          smooth=Smooth.Bezier,
          fillColor={90,34,117},
          fillPattern=FillPattern.Solid)}), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
    <p>2023, Marek Mateják</p>
</html>"));
end DropOfCommons;
