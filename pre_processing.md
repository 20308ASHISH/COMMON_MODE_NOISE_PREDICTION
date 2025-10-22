Three datasets were generated: [event, cell, ADC]

ADC1: Pedestal with mean 120

ADC2: Intrinsic noise (mean = 0, σ = 4)

ADC3: Common-mode noise (mean = 0, σ = 2) +  Intrinsic noise (mean = 0, σ = 4)

```python
    def generate_adc_data(n_events=10000, n_cells=192, pedestal_mean=120, intrinsic_sigma=4, common_sigma=2, mode='pedestal'):

        adc_data = []

        for event in range(n_events):
            if mode == 'combined':
                y = random.gauss(0, common_sigma)
            else:
                y = 0  # no common mode noise

            for cell in range(n_cells):
                if mode == 'pedestal':
                    adc = pedestal_mean
                elif mode == 'intrinsic':
                    x = random.gauss(0, intrinsic_sigma)
                    adc = pedestal_mean + x
                elif mode == 'combined':
                    x = random.gauss(0, intrinsic_sigma)
                    sign = random.choice([-1, 1])
                    adc = pedestal_mean + x + y
                else:
                    raise ValueError("mode must be 'pedestal', 'intrinsic', or 'combined'")

                adc_data.append([event, cell, adc])

        return adc_data

```

From the data obtained i saved 10 files containing ADC1 + ADC2 + ADC3.

We sorted the data based on cell...
```python
    # we take adc data and select data based on cell with events
    def compute_cell_(adc_data, n_cells=192):

        cell_values = [[] for _ in range(n_cells)]

        for event, cell, adc in adc_data:
            cell_values[cell].append(adc)

        return cell_values
```
and compute the cell data 
```python 

    def compute_cell_stats(adc_data, n_cells=192):

        cell_values = compute_cell_(adc_data, n_cells)
        # Compute mean and std for each cell
        mean_per_cell = [np.mean(vals) if vals else np.nan for vals in cell_values]
        std_per_cell  = [np.std(vals)  if vals else np.nan for vals in cell_values]

        return mean_per_cell, std_per_cell
```


and also sorted data based on events ...

```python
    # we take adc data and select data based on cell with events
    def compute_event_(adc_data, N_events=10000):

        event_values = [[] for _ in range(N_events)]

        for event, cell, adc in adc_data:
            event_values[event].append(adc)

        return event_values
```
```python

    def compute_event_stats(adc_data, N_events=192):

        event_values = compute_event_(adc_data, N_events)
        # Compute mean and std for each cell
        mean_per_event = [np.mean(vals) if vals else np.nan for vals in event_values]
        std_per_event  = [np.std(vals)  if vals else np.nan for vals in event_values]

        return mean_per_event, std_per_event
```

Then completed plots, average and standard deviation for all ADC's

After that we computed the pedestal substraction portion 
$$
ADC_{i,cell} - <ADC>_{i,cell}
$$
per channel ...
```python
    def pedestal_substraction(ADC, ADC_avg):  # give channel response only
        ADC_subtracted = []
        for  event, cell, adc in ADC:
            adc_corr = adc - ADC_avg[cell]
            ADC_subtracted.append([event, cell, adc_corr])
        return ADC_subtracted
```
we found the the ADC distribution is centered at origin with same deviation

then we compute the common mode(CM) substraction 
$$
CM(i) = \frac{1}{N_{\text{valid}}} \sum_{\text{cells}} ADC_{i,\text{cell}}
$$

The selection threshold is $3 \times \sigma_{\text{total}}^{(i)}$ for the $i^{\text{th}}$ channel.

$$
ADC_{i,cell} - <ADC>_{i,cell} - CM(i)
$$
```python
    def commom_mode_substraction(ADC, ADC_avg, ADC_deviation, N_cells=192, N_events=10000):

        # Step 1: subtract pedestal
        ADC_sub = pedestal_substraction(ADC, ADC_avg)

        # Convert to 2D array [event][cell]
        ADC_matrix = np.zeros((N_events, N_cells))
        for event, cell, adc in ADC_sub:
            ADC_matrix[event, cell] = adc

        # Step 2: compute common mode per event
        CM = np.zeros(N_events)
        for evt in range(N_events):
            valid_cells = np.abs(ADC_matrix[evt, :]) < (3 * np.array(ADC_deviation))
            CM[evt] = np.mean(ADC_matrix[evt, valid_cells])

        ADC_final = ADC_matrix - CM[:, np.newaxis]

        return ADC_final
```

WE OBSERVE THE FOLLOWING:---

### Observational Summary Table

### 📊 Observational Summary Table

| Dataset | Parameter Type | Min | Max | Mean (μ) | Std Dev (σ) |
|:--------:|:----------------------|:---:|:---:|:-----------:|:-------------:|
| **ADC2** | Mean | 119.89 | 120.094 | `ADC2_AVG` | — |
| **ADC2** | Std Dev | 3.9198 | 4.0587 | — | `ADC2_SD` |
| **ADC3** | Mean | 119.87 | 120.125 | `ADC3_AVG` | — |
| **ADC3** | Std Dev | 4.3873 | 4.5694 | — | `ADC3_SD` |
| **ADC2 (Pedestal Subtracted)** | Mean | 0.00 | 0.00 | `ADC2_AVG_PED` | — |
| **ADC2 (Pedestal Subtracted)** | Std Dev | 3.9198 | 4.0587 | — | `ADC2_SD_PED` |
| **ADC3 (Pedestal Subtracted)** | Mean | 0.00 | 0.00 | `ADC3_AVG_PED` | — |
| **ADC3 (Pedestal Subtracted)** | Std Dev | 4.3873 | 4.5694 | — | `ADC3_SD_PED` |
| **ADC2 (Common Mode Subtracted)** | Mean | 0.00 | 0.00 | `ADC2_AVG_CM` | — |
| **ADC2 (Common Mode Subtracted)** | Std Dev | 3.9138 | 4.0484 | — | `ADC2_SD_CM` |
| **ADC3 (Common Mode Subtracted)** | Mean | 0.00 | 0.00 | `ADC3_AVG_CM` | — |
| **ADC3 (Common Mode Subtracted)** | Std Dev | 3.9161 | 4.0616 | — | `ADC3_SD_CM` |
