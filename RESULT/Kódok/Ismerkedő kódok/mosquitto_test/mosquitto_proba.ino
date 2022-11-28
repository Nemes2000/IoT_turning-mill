#include <WiFi.h>               //built in library for Wifi connection
#include <Wire.h>               //built in lib for iic communication
#include <PubSubClient.h>       //for MQTT connection (https://github.com/knolleary/pubsubclient)

//###################################################|SETTINGS|##################################################################
//--------------------------------------|CAN BE SETTED|--------------------------------------------------------------------------
#define MIC_PIN 36                      //SET: the pin where you get the analog input from microphon
#define TOPIC "topic/mpu6050"           //SET: the choosen topic

IPAddress mqttServer(192,168,43,74);     //SET: to the mqtt server ip address
const char* SSID = "Nemes";              //SET: to the wifi name
const char* PWD = "nincsen123";          //SET: to the wifi password

//--------------------------------------|SHOULD NOT BE CHANGED|-------------------------------------------------------------------
WiFiClient espClient;                    //a wifi object
PubSubClient client(espClient);          //an MQTT connection object
String msg;                              //assist to send data to MQTT server
//#################################################################################################################################

void setup() {
  Serial.begin(115200);

  //set up the wifi connection
  connectToWiFi();

  //set up the MQTT server connection
  client.setServer(mqttServer, 1883);

  Serial.println("Initialization Finished!");
}

void loop() {
  //checking for Wifi connection loss
  WiFi.mode(WIFI_STA);
  if (!client.connected())
    reconnect();

  //check for MQTT server conection loss
  client.loop();

  Serial.print(F("\nSending data"));
  msg = String(analogRead(MIC_PIN));

  //sending data
  if (! client.publish(TOPIC,msg.c_str())) {
    Serial.println("Failed");
  } else {
    Serial.println("OK!");
  }
}

//############################################|FUNCTIONS FOR CONNECT VIA WIFI TO MQTT|##################################################################

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
      client.subscribe(TOPIC);
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
  Serial.print(TOPIC);
  Serial.print(". Message: ");
  String messageTemp;

  for (int i = 0; i < length; i++) {
    Serial.print((char)message[i]);
    messageTemp += (char)message[i];
  }
  Serial.println();
}