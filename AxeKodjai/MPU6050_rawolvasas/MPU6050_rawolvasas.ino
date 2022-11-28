// Basic demo for accelerometer readings from Adafruit MPU6050

#include <Adafruit_MPU6050.h>
#include <Adafruit_Sensor.h>
#include <Wire.h>

//#include <arduinoFTT.h>

#include <Wifi.h>

IPAddress mqttServer(192, 168, 43, 74);  //mqtt server ip address
const char* SSID = "Nemes";              //wifi name
const char* PWD = "nincsen123";          //wifi password

//#define WLAN_SSID       "Nemes2000"
//#define WLAN_PASS       "nincsen123"

#define AIO_SERVER      "Proba"
#define AIO_SERVERPORT  1883
#define AIO_USERNAME    ""
#define AIO_KEY         ""

Adafruit_MPU6050 mpu;
WiFiClient client;

//Adafruit_MQTT_Client mqtt(&client, AIO_SERVER, AIO_SERVERPORT, AIO_USERNAME, AIO_KEY);
//Adafruit_MQTT_Publish rezgesCsatorna = Adafruit_MQTT_Publish(&mqtt, AIO_USERNAME "/esp/rezges");


sensors_event_t a, g, temp;
unsigned long elotte, utana;
unsigned long sampleTime = 30; //Ennyi mikorsec
float adatok[256] = {0};
int ind = 0;

void MQTT_connect();  // BUG miatt kellhet

void setup(void) {
  Serial.begin(115200);
  while (!Serial) {
    delay(10);
  }

  // WIFI csatlakozas
  Serial.print("Connecting to ");
  Serial.println(WLAN_SSID);

  WiFi.begin(WLAN_SSID, WLAN_PASS);
  while (WiFi.status() != WL_CONNECTED) {
    delay(500);
    Serial.print(".");
  }
  Serial.println();
  Serial.println("WiFi connected");
  Serial.println("IP address: "); Serial.println(WiFi.localIP());


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
  Serial.println("Initialization Finished!");
  delay(100);

  elotte = micros();
}

void loop() {

  mpu.getEvent(&a, &g, &temp);

  Serial.print(F("\nSending data"));
  Serial.print(a.acceleration.x);
  Serial.print("...");
  if (! rezgesCsatorna.publish(a.acceleration.x)) {
    Serial.println(F("Failed"));
  } else {
    Serial.println(F("OK!"));
  }

  /*if(micros() > elotte + sampleTime){
    if(ind < 256){
      adatok[ind++] = a.acceleration.x;
    }
    else{
      // betelt a tomb
      Serial.print("Tomb: {");
      for(int i = 0; i < 256; i++){
        Serial.print(String(adatok[i]) + ", ");
      }
      Serial.println("}");

      ind = 0;
      // FFT

      //float fft = arduinoFFT()
    }

    elotte = micros();
  }

  /*elotte = micros();
  mpu.getEvent(&a, &g, &temp);
  utana = micros();

  Serial.println(utana-elotte);
  Serial.print("Pitch: ");
  Serial.println(g.gyro.pitch);*/
}

void MQTT_connect() {
  int8_t ret;

  // Stop if already connected.
  if (mqtt.connected()) {
    return;
  }

  Serial.print("Connecting to MQTT... ");

  uint8_t retries = 3;
  while ((ret = mqtt.connect()) != 0) { // connect will return 0 for connected
       Serial.println(mqtt.connectErrorString(ret));
       Serial.println("Retrying MQTT connection in 5 seconds...");
       mqtt.disconnect();
       delay(5000);  // wait 5 seconds
       retries--;
       if (retries == 0) {
         // basically die and wait for WDT to reset me
         while (1);
       }
  }
  Serial.println("MQTT Connected!");
}