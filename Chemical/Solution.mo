within Chemical;
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
</html>",
      info="<html>
<h4>amountOfSolution = &sum; amountOfSubstances</h4>
<h4>mass = &sum; massOfSubstances</h4>
<h4>volume = &sum; volumeOfSubstances</h4>
<h4>freeGibbsEnergy = &sum; freeGibbsEnergiesOfSubstances</h4>
<p>To calculate the sum of extensive substance's properties is misused the Modelica \"flow\" prefix even there are not real physical flows. </p>
</html>"));
end Solution;
