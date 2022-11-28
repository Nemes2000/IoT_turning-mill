#include "EspMQTTClient.h"    // MQTT lib (https://github.com/plapointe6/EspMQTTClient)

// ###################################|MQTT KAPCSOLAT KIALAKITASA|#####################################################
// WiFi beallitasok
EspMQTTClient client(
  "",                   // BSID
  "",                   // AC Jelszo
  "",                   // MQTT Szerver IP cime
  "",                   // MQTT Szerver Felhasznaloven (Ha nem kell a szerernek felhasznalao nevet megadni akkor TOROLD KI EZT A SORT)
  "",                   // MQTT Szerver Jelszava (Ha nem kell a szerernek jelszot megadni akkor TOROLD KI EZT A SORT)
  "",        // MQTT kliens neve (ESP neve, tetszoleges)
  1883                  // MQTT szerver portja
);
// ####################################################################################################################

// ########################################|MERESI BEALLITASOK|########################################################
#define FAZIS_L1_PIN 36
#define FAZIS_L2_PIN 39
#define FAZIS_L3_PIN 34
#define ELTOLAS_PIN 35
#define MQTT_TOPIC_L1 "ESP/ARAM/L1"  // Az a csatorna/topic ahova kuldi az adatot
#define MQTT_TOPIC_L2 "ESP/ARAM/L2"  // Az a csatorna/topic ahova kuldi az adatot
#define MQTT_TOPIC_L3 "ESP/ARAM/L3"  // Az a csatorna/topic ahova kuldi az adatot

int valueL1 = 0;
int valueL2 = 0;
int valueL3 = 0;
int valueEltolas = 0;

double measuredVoltEltolas = 0.0;
double measuredAmperL1 = 0.0;
double measuredAmperL2 = 0.0;
double measuredAmperL3 = 0.0;

// ##########################################|SEGED FUGGVENYEK|########################################################
void kiirMerest_DEBUG(){
  Serial.print(measuredAmperL1);
  Serial.print("A\t");
  Serial.print(measuredAmperL2);
  Serial.print("A\t");
  Serial.print(measuredAmperL3);
  Serial.println("A");
}

double valueToVolt(int value){
  const double maxSample = 4096.0;
  const double maxVoltage = 3.3;

  return (double)value/maxSample*maxVoltage;
}

double calcAmperfromRawVolt(double rawVolt){
  rawVolt -= measuredVoltEltolas;       // Ezzel vissza toljuk a 0-ba
  double rawAmper = rawVolt / 100.0;    // U/R, ahol R = 100 Ohm
  return rawAmper * 1000.0;             // Tekercs miatti valtoszam, 1000x kisebb aramot merunk mint eredetileg ezert szorozzuk vissza
}

void onConnectionEstablished()
{
  // Amikor a kapcsolat letrejott itt opcionalisan lehet valamit csinalni
  // FONTOS:
  //  Ha nem kell/akarsz csinalni semmit akkor is meg kell hagyni a fgv-t, ures torzzsel
  //  kulonben hibat fog dobni a fordito!
}

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

  // PIN-ek beallitasa
  pinMode(FAZIS_L1_PIN, INPUT);
  pinMode(FAZIS_L2_PIN, INPUT);
  pinMode(FAZIS_L3_PIN, INPUT);
  pinMode(ELTOLAS_PIN, INPUT);

  // Init vege
  Serial.println("Initialization Finished!");
}

void loop() {

  // MQTT-hez
  client.loop();

  // Mintavetelezesek
  valueEltolas = analogRead(ELTOLAS_PIN);
  valueL1 = analogRead(FAZIS_L1_PIN);
  valueL2 = analogRead(FAZIS_L2_PIN);
  valueL3 = analogRead(FAZIS_L3_PIN);

  // Veszultseg kiszamolasa az eltolasa
  // Hamarabb kell megallapitani a calcAmperfromRawVolt fgv miatt!
  measuredVoltEltolas = valueToVolt(valueEltolas);

  // Aremerossegek kiszamolasa a fazisokra
  measuredAmperL1 = calcAmperfromRawVolt(valueToVolt(valueL1));
  measuredAmperL2 = calcAmperfromRawVolt(valueToVolt(valueL2));
  measuredAmperL3 = calcAmperfromRawVolt(valueToVolt(valueL3));

  // Adatok kuldese/kiirasa
  // kiirMerest_DEBUG();
  client.publish(MQTT_TOPIC_L1, String(measuredAmperL1));
  client.publish(MQTT_TOPIC_L2, String(measuredAmperL2));
  client.publish(MQTT_TOPIC_L3, String(measuredAmperL3));
  
  //delay(100);
}