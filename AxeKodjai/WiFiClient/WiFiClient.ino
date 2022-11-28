/*
 *  This sketch sends data via HTTP GET requests to data.sparkfun.com service.
 *
 *  You need to get streamId and privateKey at data.sparkfun.com and paste them
 *  below. Or just customize this script to talk to other HTTP servers.
 *
 */

#include <WiFi.h>
#include <ArduinoMqttClient.h>

const char* ssid     = "Telekom-E4778B";
const char* password = "5nn3cxbcgtg964dn";

MqttClient mqttClient(WiFI);

const char broker[] = "test.mosquitto.org";
int        port     = 1883;
const char topic[]  = "real_unique_topic";

void setup()
{
    Serial.begin(115200);
    delay(10);

    // We start by connecting to a WiFi network

    Serial.println();
    Serial.println();
    Serial.print("Connecting to ");
    Serial.println(ssid);

    WiFi.begin(ssid, password);

    while (WiFi.status() != WL_CONNECTED) {
        delay(500);
        Serial.print(".");
    }

    Serial.println("");
    Serial.println("WiFi connected");
    Serial.println("IP address: ");
    Serial.println(WiFi.localIP());
}

int value = 0;
byte dbNet;

void loop()
{
    dbNet = WiFi.scanNetworks();
    Serial.print("Felfedezett halozatok szama: ");
    Serial.println(dbNet);
}

