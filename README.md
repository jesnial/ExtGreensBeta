# ExtGreensBeta
Extended Green's Law - Calculating Beta <br/> <br/>

Rapidly predicting tsunami wave height up to the shore is key to determining the threat it poses. On the Mediterranean coast of France in particular, there is very little warning time due to proximity of earthquake sources and the amounts of both population and infrastructure are quite consequent. Thus, unlike on the Pacific coast, where sources are further away from the coast, this area requires fast prediction tools. Nowadays, one can have a rough estimate by running a simulation over a coarse grid in the area of interest, and by then using Green's law : <img src="https://render.githubusercontent.com/render/math?math=\eta _2= \eta _1 (\frac{h_1}{h_2})^{1/4}"> <br/> 
with <img src="https://render.githubusercontent.com/render/math?math=\eta _1, \eta _2"> wave height at reference and calculation points, and <img src="https://render.githubusercontent.com/render/math?math=h _1, h _2"> the corresponding depths. However, this law does not take into account site specificities, that is why thge extended Green's law was developed. The extended Green's law relies on an empirical parameter,<img src="https://render.githubusercontent.com/render/math?math=\beta">, which tries to capture the local bathymetry's effect on the tsunami amplification. The aim in this project is to automate the calculation of <img src="https://render.githubusercontent.com/render/math?math=\beta"> parameters at each point in any given grid, so simulations based on the extended Green's law can be used for tsunami warning. This way, one would have a much more precise maximum wave height estimate, just as fast as one would get a rough one with the current methods. A fully automated <img src="https://render.githubusercontent.com/render/math?math=\beta"> calculation will be studied. This method uses bathymetry grids and previous tsunami simulation to find the optimum <img src="https://render.githubusercontent.com/render/math?math=\beta"> at each point in the grid through a gradient descent. 
Thus, running an extended Green's law caluclation using a <img src="https://render.githubusercontent.com/render/math?math=\beta"> grid will produce quick and accurate results, without using computationnally expensive simulations.

<br/>
<br/>
© Jessie Levillain 2019
