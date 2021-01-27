# Chemical
Modelica Library of Chemical processes

[1st Price, Modelica Library Award, 11th International Modelica Conference, Sept. 21-23, 2015 in Versailles, France](http://patf-biokyb.lf1.cuni.cz/wiki/_media/wiki/user/chemical_library_-_1st_price.pdf)


## Current release


Download [Chemical 1.4.0 (2021-01-27)](../../archive/v1.4.0.zip)


## Main references
 * Documentations
   * [Theoretical background](https://github.com/MarekMatejak/Chemical/raw/master/Chemical/Resources/Documentation/Chemical.pdf)
 * Dissertation thesis
   * [Fromalization of Integrative Physiology](https://github.com/MarekMatejak/dissertation/raw/master/thesis.pdf)
 
## Description

* Chemical Solution 
  * full thermodynamical state support: amount of substances, pressure, volume, temperature, electric potential, enthalpy, entropy
  * thermal, electric and mechanic connectors from Modelica Standard Library 3.2
* Chemical Substance 
  * ideal gas substance model 
    * with 6-parameters record as definition of gaseous chemical substance (all parameters of the substance are well described tabulated values)
  * incompressible substance model 
    * with 7-parameters record as definition of liquid or solid chemical substance (all parameters of the substance are well described tabulated values)
  * example of more than 30 fully defined chemical substances 
  * template for user substance models
* Chemical Reaction
  * any number of reactants and products
  * dissociation coefficient from free Gibbs energies (from substances definitions)
  * temperature dependences and heat flows from free enthalpies (from substances definitions)
  * new epoch making kinetics based with better fit with data
* Electro-chemical cell (batteries)
  * reduction-oxidation reactions with electron transfer to electic circuits
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

## License (BSD 3-Clause)

Copyright (c) 2008-2020, Marek Mateják, Charles University in Prague

All rights reserved. 

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. 
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Development and contribution
Release manager: [Marek Mateják](https://github.com/MarekMatejak)

Contributor: 
[Marek Mateják](https://github.com/MarekMatejak),
[Filip Ježek](https://github.com/filip-jezek), 
[tbeu](https://github.com/tbeu),
[dietmarw](https://github.com/dietmarw) 
