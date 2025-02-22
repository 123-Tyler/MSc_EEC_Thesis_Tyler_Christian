# READ ME 

## How to Speak Sparrow: A Modern Classification System of House Sparrow (*Passer domesticus*) Vocalisations

#### Student: Tyler Christian

#### Supervisors: Dr Julia Schroeder, Dr Ambre Salis

#### Course: MSc Ecology, Evoution and Conservation

-   code_and_analysis.R
    -   All code for statistical analysis and figure graphics used in this project.

-   behavioural_data.csv
    -   Unfiltered behavioural data collected during fieldwork.
        -   trial - Trial number (1-34)
        -   observation - Observation number, commutative across trials
        -   date - When the trial took place
        -   site - Where the trail took place
        -   audio_start - Time AudioMoth was turned on (hhmmss)
        -   dis_video_start - Time distant camera was turned on (hhmmss)
        -   close_video_start - Time close-up camera was turned on (hhmmss)
        -   angle_video_start - Time angled camera was turned on (hhmmss)
        -   audio_finish - Time AudioMoth was turned off (hhmmss)
        -   dis_video_finish - Time distant camera was turned off (hhmmss)
        -   close_video_finish - Time close-up camera was turned off (hhmmss)
        -   angle_video_finish - Time angled camera was turned off (hhmmss)
        -   audio_marker - Name of .wav file
        -   time - Time observation took place (hhmmss)
        -   number_male - Total male house sparrows present
        -   number_female - Total female house sparrows present
        -   total_sparrow - Total house sparrows present
        -   stationary - Binary data for stationary behaviour (zeros were left blank)
        -   aggressive - Binary data for aggresive behaviour (zeros were left blank)
        -   freeze - Binary data for freezing behaviour (zeros were left blank)
        -   mobbing - Binary data for mobbing behaviour (zeros were left blank)
        -   hide - Binary data for hiding behaviour (zeros were left blank)
        -   mating - Binary data for mating behaviour (zeros were left blank)
        -   no_sparrows_arrive - Number of house sparrows flying into the observation area
        -   no_sparrows_leave - Number of house sparrows leaving the observation area
        -   wild_starling - Number of wild starlings present
        -   wild_corvid - Number of wild corvids present (raven, rook, jackdaw magpie)
        -   wild_gull - Number of wild gulls present (mainly herring gulls)
        -   wild_raptor - Number of wild gulls present (only red kites at Ascot site)
        -   wild_pigeon - Number of wild pigeons present
        -   wild_goldfinch - Number of wild goldfinches present (relatively common at Lundy)
        -   wild_other_passerine - number of other wild passerine species present
        -   human_disturbance - Binary data of any form of human disturbance (e.g. talking, walking, cars, machinery...)
        -   first_reveal - Binary data for when the taxidermy box was first opened
        -   taxidermy_present - Type of taxidermy stimuli present (either: closed, empty, sparrow, woodpecker, starling, or sparrowhawk)
        -   open - Binary data for whether the taxidermy box was open or closed
        -   usable_stimuli - Observations marked with an "x" were filtered as usable for analysis
        -   notes - Additional observational data

-   avisoft_data.csv
    -   Extracted acoustic parameters of the clearest calls using Avisoft SASLab Lite version 5.3.2
    -   For further information, visit: <https://avisoft.com/tutorials/measuring-sound-parameters-from-the-spectrogram-automatically/>
        -   audio_marker - Name of .wav file
        -   site - Where the trail took place
        -   call - Initials for the type of call
        -   segment - Numbered segment containing extracted acoustic parameters (ranked per call)
        -   duration - Length of the segment (s)
        -   peak_freq_start - Peak frequency at start of the segment (Hz)
        -   peak_amp_start - Peak amplitude at start of the segment (m)
        -   band_width_start - Band Width at the start of the segment (Hz)
        -   quart_25_start - 25% quartile at the start of the segment
        -   quart_50_start - 50% quartile at the start of the segment
        -   quart_75_start - 75% quartile at the start of the segment
        -   entropy_start - Spectral entropy at the start of the segment
        -   peak_freq_end - Peak frequency at end of the segment (Hz)
        -   peak_amp_end - Peak amplitude at end of the segment (m)
        -   band_width_end - Band Width at the end of the segment (Hz)
        -   quart_25_end - 25% quartile at the end of the segment
        -   quart_50_end - 50% quartile at the end of the segment
        -   quart_75_end - 75% quartile at the end of the segment
        -   entropy_end - Spectral entropy at the end of the segment
        -   peak_freq_mean - Peak frequency at mean of the segment (Hz)
        -   peak_amp_mean - Peak amplitude at mean of the segment (m)
        -   band_width_mean - Band Width at the mean of the segment (Hz)
        -   quart_25_mean - 25% quartile at the mean of the segment
        -   quart_50_mean - 50% quartile at the mean of the segment
        -   quart_75_mean - 75% quartile at the mean of the segment
        -   entropy_mean - Spectral entropy at the mean of the segment

-   koe_data.csv
    -   Total segmented calls from each recording using the Koe web-based software.
        -   filename - Name of .wav file
        -   duration - Duration of the audio file (s)
        -   added - Date uploaded (dd/mm/yyyy)
        -   record_date - [UNUSED]
        -   sex - [UNUSED]
        -   species - [UNUSED]
        -   quality - Either "GOOD" (contains usable house sparrow calls), "NONE" (contains no usable house sparrow calls), or "INVALID" (invalid audio file)
        -   individual - Identification number addeded by the software
        -   track - Identification number addeded by the software
        -   sequence - Order of segments containing identified house sparrow calls
        -   note - [UNUSED]

-   Fig_3.1.1_PCA.jpeg
    -   PCA of extracted acoustic parameters from Avisoft.
    
-   Fig_3.3.1_before_after_plot.jpeg
    -   Line plot showing change in call frequencies before and after threatening stimuli.
    
-   Fig_3.4.1_dendrogram_optimised.jpeg
    -   Hirearchical cluster analysis of extracted acoustic parameters using data from Avisoft.
