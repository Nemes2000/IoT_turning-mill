#include <Adafruit_MPU6050.h>
#include <Adafruit_Sensor.h>
#include <Wire.h>

#include <arduinoFFT.h>

#define MINTAVETEL_DB 4096             // Csak kettes hatvany lehet
#define SAMPLING_FREQUENCY 512

Adafruit_MPU6050 mpu;                 // MPU modul
sensors_event_t a, g, temp;           // MPU altal mert adatok

double mertAdatok_Re[MINTAVETEL_DB];  // Szenzor altal mert adatok valos resze
double mertAdatok_Im[MINTAVETEL_DB];  // kepzetes resze
unsigned int mintavetelezesi_ido;

arduinoFFT FFT = arduinoFFT();        // FTT

unsigned long elotte;

void kiirTomb(double* t, int n, const char* szoveg){
  Serial.print(szoveg);
  Serial.print("{");

  for(int i = 0; i < n-1; ++i){
    Serial.print(String(t[i])+", ");
  }
  Serial.println(String(t[n-1])+"}");
}

void setup() {
  Serial.begin(115200);
  while (!Serial) {
    delay(10);
  }

  // MPU keresese
  if (!mpu.begin()) {
    Serial.println("Failed to find MPU6050 chip");
    while (1) {
      Serial.print(".");
      delay(200);
    }
  }
  // MPU beallitasa, config
  mpu.setAccelerometerRange(MPU6050_RANGE_2_G);
  mpu.setGyroRange(MPU6050_RANGE_250_DEG);
  mpu.setFilterBandwidth(MPU6050_BAND_21_HZ);

  mintavetelezesi_ido = round(1000000 * (1.0 / SAMPLING_FREQUENCY));

  // Inicializalasok vege
  Serial.println("Initialization Finished!");
  delay(100);

  elotte = micros();
}

int indx = 0;
long FFTElott, FFTUtan;

void loop() {
  // Adatok olvasasa
  mpu.getEvent(&a, &g, &temp);

  if(elotte + mintavetelezesi_ido < micros()){
    // Ekkor meg nem telt meg a tomb
    if(indx < MINTAVETEL_DB){
      mertAdatok_Re[indx] = a.acceleration.x;
      mertAdatok_Im[indx++] = 0.0;
    }
    // Ekkor telt meg, ki kell ertekelni FFT-vel
    else{
      //kiirTomb(mertAdatok_Re, MINTAVETEL_DB, "Elotte: ");
      FFTElott = micros();
      FFT.Windowing(mertAdatok_Re, MINTAVETEL_DB, FFT_WIN_TYP_RECTANGLE, FFT_FORWARD);
      FFT.Compute(mertAdatok_Re, mertAdatok_Im, MINTAVETEL_DB, FFT_FORWARD);
      FFT.ComplexToMagnitude(mertAdatok_Re, mertAdatok_Im, MINTAVETEL_DB);
      FFTUtan = micros();
      //kiirTomb(mertAdatok_Re, MINTAVETEL_DB, "FFT utan: ");

      Serial.println("***FFT SZAMITAS IDEJE*** " + String((FFTUtan-FFTElott)));

      // Grafikonva iras
      
      /*for(int i = 0; i < MINTAVETEL_DB; ++i){
        //Serial.print((i * 1.0 * SAMPLING_FREQUENCY) / MINTAVETEL_DB, 1);
        Serial.print(-50); Serial.print("\t");
        Serial.println(mertAdatok_Re[i]); Serial.print("\t");
        Serial.print(700); Serial.print("\t");
      }*/
      

      indx = 0;
    }

    elotte = micros();
  }
}
