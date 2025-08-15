# PPG_to_BP_Classification
This project implements a real-time photoplethysmography (PPG) to Blood pressure classification system on an STM32 microcontroller.
The system captures raw PPG signals, extracts features, runs them through a Random Forest Regressor model trained in Python and exported using emlearn, and displays the classification result (Hypotension / Normal / Hypertension) on an ILI9341 TFT display over SPI.

<img width="2609" height="983" alt="Proposed_model" src="https://github.com/user-attachments/assets/866ba6eb-2cd9-4955-8ec4-002421923893" />


### PPG Signal Acquisition
Captures 20 seconds of IR PPG data at 125 Hz sampling frequency (2500 samples).

### Feature Extraction
Processes the PPG waveform to compute 36 features which is given as input to the model to predict the BP class

### ML Model Inference
Model trained in Python (scikit-learn) and exported to C using emlearn.

### Result Display
Class label shown on a 2.4" / 2.8" ILI9341 SPI TFT display.
