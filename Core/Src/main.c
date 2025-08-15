#include "stm32f4xx_hal.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ppg_model.h"     // Your emlearn exported model
#include "ILI9341_STM32_Driver.h"
#include "ILI9341_GFX.h"

/* ======== User Config ======== */
#define FS_HZ                 125.0f      // sampling rate
#define HPF_CUTOFF_HZ         0.5f        // remove baseline drift
#define LPF_CUTOFF_HZ         8.0f        // smooth PPG
#define UART_BAUD             115200

// Cycle detection
#define MIN_PEAK_DISTANCE_S   0.30f       // 300 ms refractory
#define THRESH_INIT           50.0f       // initial threshold (ADC counts after filtering)
#define THRESH_DECAY          0.99f       // slow decay of adaptive threshold
#define THRESH_BOOST          0.20f       // fraction of new peak height to raise threshold

/* ======== Globals / Handles ======== */
ADC_HandleTypeDef hadc1;
TIM_HandleTypeDef htim3;
UART_HandleTypeDef huart2;

/* State for filters */
typedef struct {
    float prev_x;
    float prev_y;
    float alpha;
} OnePoleHPF;

typedef struct {
    float prev_y;
    float alpha;
} OnePoleLPF;

static OnePoleHPF hpf;
static OnePoleLPF lpf;

/* Cycle detection state */
typedef struct {
    float prev;           // previous sample (filtered)
    float dprev;          // previous derivative
    uint32_t n;           // sample index
    uint32_t last_peak_n; // last detected peak index
    float thresh;         // adaptive threshold
    uint8_t armed;        // armed after rising slope over threshold
} PeakDet;

static PeakDet pdet;

/* Scratch */
static volatile uint16_t latest_adc = 0;
static volatile uint8_t  sample_ready = 0;

/* ======== Utility: UART print ======== */
static void uprintf(const char *fmt, ...) {
    char buf[128];
    va_list ap;
    va_start(ap, fmt);
    int len = vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    if (len > 0) {
        HAL_UART_Transmit(&huart2, (uint8_t*)buf, (uint16_t)len, HAL_MAX_DELAY);
    }
}

/* ======== Filters ======== */
// High-pass (one pole): y[n] = alpha * (y[n-1] + x[n] - x[n-1])
static void HPF_Init(OnePoleHPF *f, float fs_hz, float fc_hz) {
    float T = 1.0f / fs_hz;
    float tau = 1.0f / (2.0f * (float)M_PI * fc_hz);
    f->alpha = tau / (tau + T);
    f->prev_x = 0.0f;
    f->prev_y = 0.0f;
}
static inline float HPF_Process(OnePoleHPF *f, float x) {
    float y = f->alpha * (f->prev_y + x - f->prev_x);
    f->prev_x = x;
    f->prev_y = y;
    return y;
}

// Low-pass (one pole): y[n] = y[n-1] + alpha * (x[n] - y[n-1])
static void LPF_Init(OnePoleLPF *f, float fs_hz, float fc_hz) {
    float T = 1.0f / fs_hz;
    float tau = 1.0f / (2.0f * (float)M_PI * fc_hz);
    f->alpha = T / (tau + T);
    f->prev_y = 0.0f;
}
static inline float LPF_Process(OnePoleLPF *f, float x) {
    float y = f->prev_y + f->alpha * (x - f->prev_y);
    f->prev_y = y;
    return y;
}

/* ======== Cycle Detection ======== */
static void PeakDet_Init(PeakDet *pd) {
    memset(pd, 0, sizeof(*pd));
    pd->thresh = THRESH_INIT;
    pd->last_peak_n = 0;
    pd->armed = 0;
}

/*
  Logic:
  - Compute derivative d = y[n] - y[n-1].
  - Arm when y > thresh and derivative > 0 (rising).
  - Detect a peak when derivative crosses from + to − (dprev > 0 && d <= 0) while armed,
    and min distance from last peak is satisfied.
  - Update adaptive threshold toward (1-THRESH_BOOST)*thresh + THRESH_BOOST*|peak|.
*/
static int PeakDet_Process(PeakDet *pd, float y, float *ibi_ms_out) {
    int found = 0;
    float d = y - pd->prev;

    // Refractory: enforce min distance
    uint32_t min_dist = (uint32_t)(MIN_PEAK_DISTANCE_S * FS_HZ);
    uint8_t refractory_ok = (pd->n - pd->last_peak_n) > min_dist;

    // Arm on rising over threshold
    if (!pd->armed && refractory_ok && y > pd->thresh && d > 0.0f) {
        pd->armed = 1;
    }

    // Peak when slope changes to negative
    if (pd->armed && (pd->dprev > 0.0f) && (d <= 0.0f) && refractory_ok) {
        uint32_t peak_n = pd->n;
        uint32_t ibi_samples = (pd->last_peak_n > 0) ? (peak_n - pd->last_peak_n) : 0;
        pd->last_peak_n = peak_n;
        pd->armed = 0;
        found = 1;

        // adaptive threshold
        float mag = fabsf(y);
        pd->thresh = THRESH_DECAY * pd->thresh + (1.0f - THRESH_DECAY) * mag;
        pd->thresh = pd->thresh * (1.0f + THRESH_BOOST);

        if (ibi_samples > 0 && ibi_ms_out) {
            *ibi_ms_out = (1000.0f * ibi_samples) / FS_HZ;
        }
    } else {
        // slow decay toward zero to avoid being stuck high
        pd->thresh *= THRESH_DECAY;
    }

    pd->dprev = d;
    pd->prev  = y;
    pd->n++;
    return found;
}

/* ======== CubeMX-like init prototypes ======== */
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_USART2_UART_Init(void);
static void MX_ADC1_Init(void);
static void MX_TIM3_Init(void);

/* ======== ADC callback ======== */
void HAL_ADC_ConvCpltCallback(ADC_HandleTypeDef* hadc) {
    if (hadc->Instance == ADC1) {
        latest_adc = HAL_ADC_GetValue(hadc);
        sample_ready = 1;
    }
}


// Redirect printf to UART
int __io_putchar(int ch) {
    HAL_UART_Transmit(&huart2, (uint8_t *)&ch, 1, HAL_MAX_DELAY);
    return ch;
}






void display_result(int result)
{
    ILI9341_Fill_Screen(WHITE);
    ILI9341_Set_Rotation(SCREEN_VERTICAL_1);
    ILI9341_Set_Text_Size(2);
    ILI9341_Set_Text_Colour(BLACK, WHITE);

    if (result == 0) {
        ILI9341_Print_Text("Hypotension", 10, 50, BLACK, WHITE, 2);
    } else if (result == 1) {
        ILI9341_Print_Text("Normal BP", 10, 50, BLACK, WHITE, 2);
    } else if (result == 2) {
        ILI9341_Print_Text("Hypertension", 10, 50, BLACK, WHITE, 2);
    } else {
        ILI9341_Print_Text("Unknown", 10, 50, BLACK, WHITE, 2);
    }
}




/* ======== Main ======== */
int main(void)
{
    HAL_Init();
    SystemClock_Config();
    MX_GPIO_Init();
    MX_USART2_UART_Init();
    MX_ADC1_Init();
    MX_TIM3_Init();

    // Initialize filters & peak detector
    HPF_Init(&hpf, FS_HZ, HPF_CUTOFF_HZ);
    LPF_Init(&lpf, FS_HZ, LPF_CUTOFF_HZ);
    PeakDet_Init(&pdet);

    // Start timer -> ADC (TRGO), ADC interrupt mode
    HAL_TIM_Base_Start(&htim3);
    HAL_ADC_Start_IT(&hadc1);

//    printf("\r\nPPG Acquisition @ %.1f Hz, HPF=%.2f Hz, LPF=%.2f Hz\r\n",
//           FS_HZ, HPF_CUTOFF_HZ, LPF_CUTOFF_HZ);

    while (1)
    {
        if (sample_ready)
        {
            sample_ready = 0;

            // Convert ADC reading to float
            float x = (float)latest_adc;

            // Filtering: HPF -> LPF
            float y = HPF_Process(&hpf, x);
            y = LPF_Process(&lpf, y);

            // Peak detection for cycle cutting
            float ibi_ms = 0.0f;
            if (PeakDet_Process(&pdet, y, &ibi_ms))
            {
                printf("Peak @ n=%lu, t=%lu ms",
                       (unsigned long)pdet.last_peak_n,
                       (unsigned long)(1000UL * pdet.last_peak_n / (uint32_t)FS_HZ));

                if (ibi_ms > 0.0f)
                {
//                    printf(", IBI=%.1f ms, HR=%.1f bpm",
//                           ibi_ms, 60000.0f / ibi_ms);
                }
                printf("\r\n");

                // Here you can add code to store the cycle samples
                // between last_peak_n and current peak for later analysis
            }
        }
        __WFI(); // Low-power wait for interrupt
    }
}






/* ======== Clock: 84 MHz (default Nucleo F401RE) ======== */
void SystemClock_Config(void)
{
    RCC_OscInitTypeDef RCC_OscInitStruct = {0};
    RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

    __HAL_RCC_PWR_CLK_ENABLE();
    __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE2);

    RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
    RCC_OscInitStruct.HSIState = RCC_HSI_ON;
    RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
    RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
    RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
    RCC_OscInitStruct.PLL.PLLM = 16;
    RCC_OscInitStruct.PLL.PLLN = 336;
    RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV4;   // 84 MHz
    RCC_OscInitStruct.PLL.PLLQ = 7;
    HAL_RCC_OscConfig(&RCC_OscInitStruct);

    RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK | RCC_CLOCKTYPE_SYSCLK |
                                  RCC_CLOCKTYPE_PCLK1 | RCC_CLOCKTYPE_PCLK2;
    RCC_ClkInitStruct.SYSCLKSource   = RCC_SYSCLKSOURCE_PLLCLK;
    RCC_ClkInitStruct.AHBCLKDivider  = RCC_SYSCLK_DIV1;
    RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;  // 42 MHz
    RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;  // 84 MHz
    HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2);
}

/* ======== ADC1 on PA0 (A0) ======== */
static void MX_ADC1_Init(void)
{
    __HAL_RCC_ADC1_CLK_ENABLE();

    ADC_ChannelConfTypeDef sConfig = {0};
    hadc1.Instance = ADC1;

    hadc1.Init.ClockPrescaler        = ADC_CLOCK_SYNC_PCLK_DIV4; // PCLK2=84MHz -> 21MHz ADC clk
    hadc1.Init.Resolution            = ADC_RESOLUTION_12B;
    hadc1.Init.ScanConvMode          = ADC_SCAN_DISABLE;
    hadc1.Init.ContinuousConvMode    = DISABLE;  // we use external trigger
    hadc1.Init.DiscontinuousConvMode = DISABLE;
    hadc1.Init.ExternalTrigConvEdge  = ADC_EXTERNALTRIGCONVEDGE_RISING;
    hadc1.Init.ExternalTrigConv      = ADC_EXTERNALTRIGCONV_T3_TRGO; // TIM3 TRGO
    hadc1.Init.DataAlign             = ADC_DATAALIGN_RIGHT;
    hadc1.Init.NbrOfConversion       = 1;
    hadc1.Init.DMAContinuousRequests = DISABLE;
    hadc1.Init.EOCSelection          = ADC_EOC_SINGLE_CONV;
    HAL_ADC_Init(&hadc1);

    sConfig.Channel      = ADC_CHANNEL_0;      // PA0
    sConfig.Rank         = 1;
    sConfig.SamplingTime = ADC_SAMPLETIME_480CYCLES; // long sample for noisy sensors
    HAL_ADC_ConfigChannel(&hadc1, &sConfig);
}

/* ======== TIM3 TRGO @ 125 Hz ======== */
static void MX_TIM3_Init(void)
{
    __HAL_RCC_TIM3_CLK_ENABLE();

    // We want update event at 125 Hz.
    // Timer clock on APB1: APB1=42 MHz, but TIMx clock is 84 MHz (x2 when APB1 prescaler != 1).
    // Compute PSC and ARR such that: f = 84e6 / ((PSC+1)*(ARR+1)) = 125
    // Choose PSC = 8399 -> timer tick = 84e6/8400 = 10000 Hz
    // Then ARR = 10000/125 - 1 = 79
    htim3.Instance = TIM3;
    htim3.Init.Prescaler         = 8399;
    htim3.Init.CounterMode       = TIM_COUNTERMODE_UP;
    htim3.Init.Period            = 79;
    htim3.Init.ClockDivision     = TIM_CLOCKDIVISION_DIV1;
    htim3.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
    HAL_TIM_Base_Init(&htim3);

    // TRGO on update
    TIM_MasterConfigTypeDef sMasterConfig = {0};
    sMasterConfig.MasterOutputTrigger = TIM_TRGO_UPDATE;
    sMasterConfig.MasterSlaveMode     = TIM_MASTERSLAVEMODE_DISABLE;
    HAL_TIMEx_MasterConfigSynchronization(&htim3, &sMasterConfig);
}

/* ======== USART2 on ST-Link VCP ======== */
static void MX_USART2_UART_Init(void)
{
    __HAL_RCC_USART2_CLK_ENABLE();

    huart2.Instance          = USART2;
    huart2.Init.BaudRate     = UART_BAUD;
    huart2.Init.WordLength   = UART_WORDLENGTH_8B;
    huart2.Init.StopBits     = UART_STOPBITS_1;
    huart2.Init.Parity       = UART_PARITY_NONE;
    huart2.Init.Mode         = UART_MODE_TX_RX;
    huart2.Init.HwFlowCtl    = UART_HWCONTROL_NONE;
    huart2.Init.OverSampling = UART_OVERSAMPLING_16;
    HAL_UART_Init(&huart2);
}

/* ======== GPIO (PA0 analog, USART pins set by Cube) ======== */
static void MX_GPIO_Init(void)
{
    __HAL_RCC_GPIOA_CLK_ENABLE();
    GPIO_InitTypeDef GPIO_InitStruct = {0};

    // PA0 as analog
    GPIO_InitStruct.Pin  = GPIO_PIN_0;
    GPIO_InitStruct.Mode = GPIO_MODE_ANALOG;
    GPIO_InitStruct.Pull = GPIO_NOPULL;
    HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);
}

/* ======== IRQ Handlers (if needed depending on HAL config) ======== */
void SysTick_Handler(void) {
    HAL_IncTick();
}
