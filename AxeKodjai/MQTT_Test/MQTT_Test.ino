#include "EspMQTTClient.h"

#include <Adafruit_MPU6050.h>
#include <Adafruit_Sensor.h>
#include <Wire.h>

// Otthoni beallitasok
/*EspMQTTClient client(
  "Telekom-E4778B",
  "5nn3cxbcgtg964dn",
  "192.168.1.228",  // MQTT Broker server ip
  "ESP32-Rezgesmero",     // Client name that uniquely identify your device
  1883              // The MQTT port, default to 1883. this line can be omitted
);*/

// BME I4.0 beallitasok
EspMQTTClient client(
  "I40TK-student",
  "I40Hallgat22",
  "172.18.0.1",  // MQTT Broker server ip
  "ESP32-Rezgesmero",     // Client name that uniquely identify your device
  1883              // The MQTT port, default to 1883. this line can be omitted
);

Adafruit_MPU6050 mpu;
sensors_event_t a, g, temp;

void setup()
{
  Serial.begin(115200);
  while (!Serial) {
    delay(10);
  }

  // MQTT beallitasok (opcionalisak)
  client.enableDebuggingMessages();
  client.enableHTTPWebUpdater();
  client.enableOTA();
  client.enableLastWillMessage("TestClient/lastwill", "I am going offline");

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

  // Init vege
  Serial.println("Initialization Finished!");
  delay(100);
}

void onConnectionEstablished()
{
  client.subscribe("esp/rezges", [](const String & payload) {
    Serial.println(payload);
  });

  client.publish("esp/chat", "This is a message");
}

void loop()
{
  mpu.getEvent(&a, &g, &temp);

  client.loop();
  client.publish("esp/rezges/x", String(a.acceleration.x));
  client.publish("esp/rezges/y", String(a.acceleration.y));
  client.publish("esp/rezges/z", String(a.acceleration.z));
  
  delay(100);
}