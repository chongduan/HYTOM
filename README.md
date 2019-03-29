

# HYTOM
Hybrid T One and Magnetization transfer (HYTOM) is a gadolinium-free cardiac magnetic resonance imaging method for detecting and characterizing myocardial infarction. Further details are described in <a href="https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27636">this paper</a>.

This repository contains all simulation codes used in the above paper. These simulations are built upon the <a href="https://github.com/mriphysics/EPG-X">EPG-X framework</a> developed by <a href="https://github.com/shaihanmalik">Shaihan Malik, King’s College London</a>.

Feel free to reach me at <a href="mailto:duanchong520@gmail.com">duanchong520@gmail.com</a> if you have any questions or want to contribute to these simulation code.

If you find this helpful, please cite <a href="https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27636">Duan et al. MRM 2018</a>.

## Requirements
<li>Matlab R2017a (earlier version may also work, not tested yet)</li>
<li>Parallel Computing Toolbox</li>

## Example simulations
### test_figure_2a_2b
2a: Simulate the bSSFP signal evolution for three scenarios: 1) single pool; 2) two pools; 3) two pools with magnetization transfer preparation.

<img src="https://raw.githubusercontent.com/chongduan/HYTOM/master/Images/bSSFP_profile.png" width="500">

2b: Plot the off-resonance profiles for each of the scenarios at the center of kspace (you may choose any echoes).

<img src="https://raw.githubusercontent.com/chongduan/HYTOM/master/Images/bSSFP_to_ss.png" width="500">


### test_figure_2c
2c: Simulate the inversion recovery curve for MOLLI and HYTOM, and fit the simulated signal to a 3-parameter model

<img src="https://raw.githubusercontent.com/chongduan/HYTOM/master/Images/Relax.png" width="500">

## License

<a href="https://choosealicense.com/licenses/mit/">MIT</a>
