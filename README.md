# DaphniaMagnaTracking

The following script is used to analyse the behaviour of Daphnia from videos taken in the laboratory after being exposed to different concentrations of a pesticide.

The original videos have undergone two operations :

- A mask was applied. 
- Reduction of the size to 25%. 


This script uses mainly three packages :

“TrackDem” (TD) : motion analysis (https://github.com/marjoleinbruijning/trackdem).

trajr (TR) : Stopping time and sinuosity analysis (McLean DJ, Skowron Volponi MA. trajr: An R package for characterisation of animal trajectories. Ethology. 2018;00:19. https://doi.org/10.1111/eth.12739.).

av (AV) : cutting videos into images (https://CRAN.R-project.org/package=av).

## Target :

The behavioural data obtained is analysed statistically to see the impact of the pesticide on the behaviour of Daphnia.

## Example of results : 

### Final image :
<img src="https://user-images.githubusercontent.com/100361637/213686593-5da1a400-bc41-4835-a9f4-d62860f3d563.png" height="500">

### Animation :

![ezgif com-gif-maker](https://user-images.githubusercontent.com/100361637/213685524-7df0e820-b234-4953-b7bd-5609b568199d.gif)
