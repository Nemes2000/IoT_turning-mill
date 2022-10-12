#include <WiFi.h>
#include <PubSubClient.h>
#include <Wire.h>

IPAddress mqttServer(192, 168, 43, 74);  //mqtt server ip address
const char* SSID = "Nemes";              //wifi name
const char* PWD = "nincsen123";          //wifi password

WiFiClient espClient;
PubSubClient client(espClient);

void setup() {
  Serial.begin(9600);

  connectToWiFi();
  client.setServer(mqttServer, 1883);
}

void loop() {
  WiFi.mode(WIFI_STA);
  if (!client.connected())
    reconnect();

  client.loop();

  long now = millis();
  long lastime = 0;
  if (now - lastime > 5000) {
    Serial.println("Publishing data..");
    client.publish("topic/proba", "proba");
    lastime = now;
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