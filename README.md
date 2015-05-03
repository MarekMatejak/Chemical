# Chemical
Modelica Library of Chemical processes

* Chemical Solution 
  * full thermodynamical state support: amount of substances, pressure, volume, temperature, electric potential, enthalpy, entropy
  * thermal, electric and mechanic connectors from Modelica Standard Library 3.2
* Chemical Substance 
  * ideal gas substance model with 5-parameter record as definition of gaseous chemical substance
  * incompressible substance model with 7-parameter record as definitions of liquid or solid chemical substance
  * example of more than 30 fully defined chemical substances 
  * template for user substance models
* Chemical Reaction
  * any number of reactants and products
  * dissociation coefficient from free Gibbs energies (from substances definitions)
  * temperature dependences and heat flows from free enthalpies (from substances definitions)
  * new epoch making kinetics based with better fit with data
* Electro-chemical cell (batteries)
  * reduction-oxidation reactions with redirection of electron flows to electic circuits
  * electrodes as solid chemical solutions
  * realistic disscharging curve without any lookup data
* Gas Dissolution
  * Henry's, Raoult's and Sieverts' in one component
  * dissolution coefficients from definition of substances
* Membrane
  * semipermeable membrane with selective transport of each substance 
  * different substances can cross the membrane in both directions at the same time
  * osmotic transport for uncharges substances to reach osmotic equilibrium
  * reach Donnan's equilibrium of electorlytes
  * generate Nernst (Goldman-Hodgkin) electric potential on membrane
* Chemical Speciation of macromolecules
  * can rapidly simplify the equilibrium on macromolecule
  * allosteric effects
  * example of hemoglobin oxygen binding
* Diffusion
* Substance flow in stream of the solution
* Degradation
* Clearance
* ...

* Based on equilibration of the electro-chemical potentials. 
* With full support of thermodynamics of chemical solutions.

