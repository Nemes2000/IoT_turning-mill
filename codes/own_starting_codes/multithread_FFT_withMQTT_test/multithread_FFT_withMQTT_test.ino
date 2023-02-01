#include <WiFi.h>
#include <PubSubClient.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_MPU6050.h>
#include <Wire.h>
#include <arduinoFFT.h>

#define MINTAVETEL_DB 4096        //nincs eleg hely -> hang kimegy lehet tobb     
#define SAMPLING_FREQUENCY 512

IPAddress mqttServer(192, 168, 43, 74);  //mqtt server ip address
const char* SSID = "Nemes";              //wifi name
const char* PWD = "nincsen123";          //wifi password

Adafruit_MPU6050 mpu;
WiFiClient espClient;
PubSubClient client(espClient);

sensors_event_t a, g, temp;
unsigned long elotte;
double mertAdatok_Re1[MINTAVETEL_DB];  // Szenzor altal mert adatok valos resze
double mertAdatok_Im1[MINTAVETEL_DB];  // kepzetes resze
unsigned int mintavetelezesi_ido;
TaskHandle_t Task1;

arduinoFFT FFT = arduinoFFT();        // FTT

void setup() {
  Serial.begin(115200);

  connectToWiFi();
  client.setServer(mqttServer, 1883);

  // MPU keresese
  if (!mpu.begin()) {
    Serial.println("Failed to find MPU6050 chip");
    while (1) {
      delay(10);
    }
  }

  // MPU beallitasa, config
  mpu.setAccelerometerRange(MPU6050_RANGE_2_G);
  mpu.setGyroRange(MPU6050_RANGE_250_DEG);
  mpu.setFilterBandwidth(MPU6050_BAND_21_HZ);

  mintavetelezesi_ido = round(1000000 * (1.0 / SAMPLING_FREQUENCY));

  Serial.println("Initialization Finished!");
  delay(100);

  elotte = micros();
}

long FFTElott, FFTUtan;
char msg_out[20];   //float publikalasahoz

void loop() {
  WiFi.mode(WIFI_STA);
  if (!client.connected())
    reconnect();

  client.loop();

  mpu.getEvent(&a, &g, &temp);

  if(elotte + mintavetelezesi_ido < micros()){
    if(indx < MINTAVETEL_DB){
      mertAdatok_Re1[indx] = a.acceleration.x;
      mertAdatok_Im1[indx++] = 0.0;
    }
    else{
       //kiirTomb(mertAdatok_Re, MINTAVETEL_DB, "Elotte: ");
      FFTElott = micros();
      FFT.Windowing(mertAdatok_Re1, MINTAVETEL_DB, FFT_WIN_TYP_RECTANGLE, FFT_FORWARD);
      FFT.Compute(mertAdatok_Re1, mertAdatok_Im1, MINTAVETEL_DB, FFT_FORWARD);
      FFT.ComplexToMagnitude(mertAdatok_Re1, mertAdatok_Im1, MINTAVETEL_DB);
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
    }
  }

  Serial.print(F("\nSending data"));
  Serial.print(a.acceleration.x);
  Serial.print("...");
  if (! client.publish("topic/mpu6050",dtostrf(a.acceleration.x,10,10,msg_out))) {
    Serial.println("Failed");
  } else {
    Serial.println("OK!");
  }
}

void connectToWiFi() {
  Serial.print("Connectiog to ");

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
      client.subscribe("topic/proba");
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