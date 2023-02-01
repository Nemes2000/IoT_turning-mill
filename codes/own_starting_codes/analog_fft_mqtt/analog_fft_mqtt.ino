#include <WiFi.h>
#include <PubSubClient.h>
#include <arduinoFFT.h>
#include <Wire.h>

#define SAMPLE_QUANTITY 4096  
#define RESULT_QUANTITY 20
#define MESSAGE_LENGTH 20

const int microphonPin = 36;

IPAddress mqttServer(192,168,43,74);    //mqtt server ip address
const char* SSID = "Nemes";             //wifi name
const char* PWD = "nincsen123";         //wifi password
WiFiClient espClient;
PubSubClient client(espClient);

double samples_Re[SAMPLE_QUANTITY];     // collected results real part by the sensor 
double samples_Im[SAMPLE_QUANTITY];     // imaginary part
int nextSapleIndex = 0;
arduinoFFT FFT = arduinoFFT();          // FTT
unsigned long beforeCollection;
unsigned long beforeFFT;
char msg_out[MESSAGE_LENGTH];           //for data publishing

void setup() {
  Serial.begin(115200);

  while(!Serial);
   Serial.println("Ready");

  connectToWiFi();
  client.setServer(mqttServer, 1883);
  beforeCollection = micros();
}

void loop() {
  WiFi.mode(WIFI_STA);
  if (!client.connected())
    reconnect();

  client.loop();

  fft();
}

void fft(){
  
  if(nextSapleIndex < SAMPLE_QUANTITY){
    samples_Re[nextSapleIndex] = analogRead(microphonPin);
    samples_Im[nextSapleIndex++] = 0.0;
  }
  else{
    unsigned long CollectingTime = micros()-beforeCollection;
    double CollectingTimeD = CollectingTime;

    Serial.println("***COLLECTING TIME*** " + String(CollectingTime));

    beforeFFT = micros();

    //FFT calculation
    FFT.Windowing(samples_Re, SAMPLE_QUANTITY, FFT_WIN_TYP_RECTANGLE, FFT_FORWARD);
    FFT.Compute(samples_Re, samples_Im, SAMPLE_QUANTITY, FFT_FORWARD);
    FFT.ComplexToMagnitude(samples_Re, samples_Im, SAMPLE_QUANTITY);
    unsigned long FFTime = micros()-beforeFFT;
    Serial.println("***FFT ALGORTIHM TIME*** " + String(FFTime));

    //sampling frequency calc
    double CollectTimeInSec = CollectingTimeD / 1000000.0;
    int sampling_feq = 4096.0 / CollectTimeInSec; // => around 5300

    // collect the highest peaks
    int results[RESULT_QUANTITY];
    MajorPeaks(samples_Re, SAMPLE_QUANTITY, sampling_feq, results); //5300 =>  4096 / collecting futas ideje sec-ben

    //sending peaks for mqtt server
    for(int i = 0; i < RESULT_QUANTITY; i++){
      Serial.print("Freq: "); 
      Serial.print(results[i]); 

      if (! client.publish("topic/microphon", dtostrf(results[i], MESSAGE_LENGTH - 3, 3, msg_out))) {
        Serial.println("Failed");
      } else {
        Serial.println("OK!");
      }
    }

    nextSapleIndex = 0;
    beforeCollection = micros();
  }
}

void MajorPeaks(double *vD, uint16_t samples, double samplingFrequency, int *peaks)
{
  int peaksYIndex[RESULT_QUANTITY];
  for(int k = 0; k < RESULT_QUANTITY; k++){
    peaksYIndex[k] = -1;
    peaks[k] = 0;
  }

	for(int ind = 0; ind < RESULT_QUANTITY; ind++){
    double maxY = 0;
    uint16_t IndexOfMaxY = 0;
    //If sampling_frequency = 2 * max_frequency in signal,
    //value would be stored at position samples/2
    for (uint16_t i = 1; i < ((samples >> 1) + 1); i++) {
      if ((vD[i-1] < vD[i]) && (vD[i] > vD[i+1])) {

        bool has_been = false;
        for(int pi = 0; pi < RESULT_QUANTITY; pi++)
          if(i == peaksYIndex[pi])
            has_been = true;

        if (vD[i] > maxY && !has_been) {
          maxY = vD[i];
          IndexOfMaxY = i;
        }
      }
    }
    double delta = 0.5 * ((vD[IndexOfMaxY-1] - vD[IndexOfMaxY+1]) / (vD[IndexOfMaxY-1] - (2.0 * vD[IndexOfMaxY]) + vD[IndexOfMaxY+1]));
    double interpolatedX = ((IndexOfMaxY + delta)  * samplingFrequency) / (samples-1);
    if(IndexOfMaxY==(samples >> 1)) //To improve calculation on edge values
      interpolatedX = ((IndexOfMaxY + delta)  * samplingFrequency) / (samples);
    // returned value: interpolated frequency peak apex
    peaksYIndex[ind] = IndexOfMaxY;
    peaks[ind] = interpolatedX;
  }
}

void connectToWiFi() {
  Serial.print("Connecting to ");

  WiFi.begin(SSID, PWD);
  Serial.println(SSID);
  while (WiFi.status() != WL_CONNECTED) {
    Serial.print(".");
    delay(500);
  }
  Serial.print("Connected.");

  Serial.println("IP address: ");
  Serial.println(WiFi.localIP());
}

void reconnect() {
  // Loop until we're reconnected
  while (!client.connected()) {
    Serial.print("Attempting MQTT connection...");
    
    // Attempt to connect
    if (client.connect("ESP32Client")) {
      Serial.println("connected");
      // Subscribe
      client.subscribe("topic/mpu9250");
    } else {
      Serial.print("failed, rc=");
      Serial.print(client.state());
      Serial.println(" try again in 5 seconds");
      Serial.println(WiFi.localIP());
      // Wait 5 seconds before retrying
      delay(5000);
    }
  }
}

//giving orders via mqtt -> do not need now
void callback(char* topic, byte* message, unsigned int length) {
  Serial.print("Message arrived on topic: ");
  Serial.print(topic);
  Serial.print(". Message: ");
  String messageTemp;

  for (int i = 0; i < length; i++) {
    Serial.print((char)message[i]);
    messageTemp += (char)message[i];
  }
  Serial.println();
}
