#include "EspMQTTClient.h"      // MQTT lib (https://github.com/plapointe6/EspMQTTClient)
#include <Wire.h>
#include <Adafruit_MLX90614.h>  // Segedlet hozza (https://create.arduino.cc/projecthub/pibots555/how-to-connect-dht11-sensor-with-arduino-uno-f4d239)

// ###################################|MQTT KAPCSOLAT KIALAKITASA|#####################################################
// WiFi beallitasok
EspMQTTClient client(
  "",                   // BSID
  "",                   // AC Jelszo
  "",                   // MQTT Szerver IP cime
  "",                   // MQTT Szerver Felhasznaloven (Ha nem kell a szerernek felhasznalao nevet megadni akkor TOROLD KI EZT A SORT)
  "",                   // MQTT Szerver Jelszava (Ha nem kell a szerernek jelszot megadni akkor TOROLD KI EZT A SORT)
  "ESP32-Homerseklet",        // MQTT kliens neve (ESP neve, tetszoleges)
  1883                  // MQTT szerver portja
);
// ########################################|MERESI BEALLITASOK|########################################################
#define MQTT_TOPIC_C "ESP/Homerseklet/Celsius"                    // Az a csatorna/topic ahova kuldi az adatot
#define MQTT_TOPIC_F "ESP/Homerseklet/Fahrenheit"                  // Az a csatorna/topic ahova kuldi az adatot

Adafruit_MLX90614 mlx = Adafruit_MLX90614();                      // MLX szenzor objektuma

double homerseklet_C, homerseklet_F;

// ########################################|SETUP|#####################################################################
void setup() {
  // Serial bekapcsolasa, megvarasa
  Serial.begin(115200);
  while (!Serial) {
    delay(10);
  }

  // MQTT beallitasok (opcionalisak)
  client.enableDebuggingMessages();
  client.enableHTTPWebUpdater();
  client.enableOTA();
  client.setMqttReconnectionAttemptDelay(1500);

  // MLX elinditasa
  mlx.begin();

  // Init vege
  Serial.println("Initialization Finished!");
}

void onConnectionEstablished()
{
  // Amikor a kapcsolat letrejott itt opcionalisan lehet valamit csinalni
  // FONTOS:
  //  Ha nem kell/akarsz csinalni semmit akkor is meg kell hagyni a fgv-t, ures torzzsel
  //  kulonben hibat fog dobni a fordito!
}

// ########################################|LOOP|######################################################################
void loop() {
  // MQTT-hez
  client.loop();

  homerseklet_C = mlx.readObjectTempC();
  homerseklet_F = mlx.readObjectTempF();

  // Adatok kuldese/kiirasa
  Serial.print("Object temp (C): ");  Serial.println(homerseklet_C);
  Serial.print("Object temp (F): ");  Serial.println(homerseklet_F);
  client.publish(MQTT_TOPIC_C, String(homerseklet_C));
  client.publish(MQTT_TOPIC_F, String(homerseklet_F));

  delay(100);
}
