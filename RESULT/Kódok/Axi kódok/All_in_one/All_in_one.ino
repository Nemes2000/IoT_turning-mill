#include "EspMQTTClient.h"    // MQTT lib (https://github.com/plapointe6/EspMQTTClient)
#include <Wire.h>
#include <Adafruit_MLX90614.h>  // Segedlet hozza (https://create.arduino.cc/projecthub/pibots555/how-to-connect-dht11-sensor-with-arduino-uno-f4d239)

// ###################################|MQTT KAPCSOLAT KIALAKITASA|#####################################################
// WiFi beallitasok
EspMQTTClient client(
  "I40TK-office",                   // BSID
  "Ipar4Ir0d4",                     // AC Jelszo
  "172.22.0.79",                   // MQTT Szerver IP cime
  "Esp-All_in_one",        // MQTT kliens neve (ESP neve, tetszoleges)
  1883                     // MQTT szerver portja
);
// ####################################################################################################################

// ########################################|MERESI BEALLITASOK|########################################################
#define FAZIS_L1_PIN 33 //sarga
#define FAZIS_L2_PIN 35 //kek
#define FAZIS_L3_PIN 34 //zold
#define ELTOLAS_PIN 18
#define MQTT_TOPIC_L1 "topic/ARAM-L1"  // Az a csatorna/topic ahova kuldi az adatot
#define MQTT_TOPIC_L2 "topic/ARAM-L2"  // Az a csatorna/topic ahova kuldi az adatot
#define MQTT_TOPIC_L3 "topic/ARAM-L3"  // Az a csatorna/topic ahova kuldi az adatot
#define MQTT_TOPIC_C "topic/Homerseklet-Celsius"                    // Az a csatorna/topic ahova kuldi az adatot
#define MQTT_TOPIC_F "topic/Homerseklet-Fahrenheit"                  // Az a csatorna/topic ahova kuldi az adatot
#define OPTO_PIN 21           // Opto pinje
#define TIMEOUT 0.4           // Hany masodperc semmitteves utan timeout-oljon
#define FINOM_EJTES true      // Ha timeout-ol a fordulatszam mero akkor fokozatosan eljtse-e a jelet vagy sem
#define FINOM_SZORZO 0.4      // Ha finom ejtest alkalmazunk, akkor milyen gyorsan csengjen le a jel
#define MQTT_TOPIC "topic/RPM"  // Az a csatorna/topic ahova kuldi az adatot

long elotte, most;            // Az elozo es mostani interrupt ideje microsec-ben
int tombIndex = 0;            // a lyukakKozottiIdo indexe, ora modszerrel (opre) toltjuk fel/frissitjuk az ertekeket a tombben
long tombOsszeg = 0.0;        // Ebbe szamoljuk mindig a tombnek az osszeget
double RPM = 0.0;             // Fordulatszam
bool timeoutVan = false;      // Van-e timeout?
long lyukakKozottiIdo[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};    // Ablakkos modszer miatt, eltaroljuk a lyukak kozotti erteket es az alapjan szamolunk
 

int valueL1 = 0;
int valueL2 = 0;
int valueL3 = 0;
int valueEltolas = 0;

Adafruit_MLX90614 mlx = Adafruit_MLX90614();                      // MLX szenzor objektuma

double homerseklet_C, homerseklet_F;
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

// Interupt fgv, minden egyes opto olvassakor meghivodik
void IRAM_ATTR interruptFgv(){
  most = micros();

  // Megfelelo helyre beirjuk a ket lyuk kozotti erteket
  if(tombIndex > 5){
    tombIndex = 0;
  }
  lyukakKozottiIdo[tombIndex++] = most-elotte;

  // Korulfordulasi idot kiszamoljuk, majd az RPM-et kiszamoljuk
  tombOsszeg = lyukakKozottiIdo[0] + lyukakKozottiIdo[1] + lyukakKozottiIdo[2] + lyukakKozottiIdo[3] + lyukakKozottiIdo[4] + lyukakKozottiIdo[5];
  RPM = 60000000.0 / (float)tombOsszeg;

  elotte = most;
}

// ####################################################################################################################

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

  // Az opto kapu pinjet inputra rakjuk
  pinMode(OPTO_PIN, INPUT);
  // Interrupt beallitasa
  attachInterrupt(digitalPinToInterrupt(OPTO_PIN), interruptFgv, RISING);

  // MLX elinditasa
  Wire.begin(23, 18);
  mlx.begin(90, &Wire);

  // PIN-ek beallitasa
  pinMode(FAZIS_L1_PIN, INPUT);
  pinMode(FAZIS_L2_PIN, INPUT);
  pinMode(FAZIS_L3_PIN, INPUT);
  pinMode(ELTOLAS_PIN, INPUT);
  
  elotte = micros();

  // Init vege
  Serial.println("Initialization Finished!");
}

void loop() {
  // Mintavetelezesek
  valueEltolas = analogRead(ELTOLAS_PIN);
  valueL1 = analogRead(FAZIS_L1_PIN);
  valueL2 = analogRead(FAZIS_L2_PIN);
  valueL3 = analogRead(FAZIS_L3_PIN);
  homerseklet_C = mlx.readObjectTempC();
  homerseklet_F = mlx.readObjectTempF();

  // Veszultseg kiszamolasa az eltolasa
  // Hamarabb kell megallapitani a calcAmperfromRawVolt fgv miatt!
  measuredVoltEltolas = valueToVolt(valueEltolas);

  // Aremerossegek kiszamolasa a fazisokra
  measuredAmperL1 = calcAmperfromRawVolt(valueToVolt(valueL1));
  measuredAmperL2 = calcAmperfromRawVolt(valueToVolt(valueL2));
  measuredAmperL3 = calcAmperfromRawVolt(valueToVolt(valueL3));

  // Van-e timeout
    if(micros() - elotte > TIMEOUT * 1000000.0){
      timeoutVan = true;
      // Lenullazzuk a tombot
      for(int i = 0; i < 6; ++i){
        lyukakKozottiIdo[i] = 0.0;
      }
    }
    else{
      timeoutVan = false;
    }

    // Ha timeout van akkor nullazni kell az RPM-et Timeout => Kikapcsoltak az esztergat
    if(timeoutVan && RPM > 0.0){

      if(FINOM_EJTES){
        if(RPM > 10.0)
          RPM -= RPM * FINOM_SZORZO;
        else
          RPM = 0.0;      
      }
      else{
        RPM = 0.0;
      }
    }

    // Adatok kuldese/kiirasa
  Serial.println(RPM);
  client.loop();
  client.publish(MQTT_TOPIC, String(RPM));

  // Adatok kuldese/kiirasa
  // kiirMerest_DEBUG();
  Serial.print("Object temp (C): ");  Serial.println(homerseklet_C);
  Serial.print("Object temp (F): ");  Serial.println(homerseklet_F);
  client.loop();
  client.publish(MQTT_TOPIC_C, String(homerseklet_C));
  client.loop();
  client.publish(MQTT_TOPIC_F, String(homerseklet_F));
  client.loop();
  client.publish(MQTT_TOPIC_L1, String(measuredAmperL1));
  client.loop();
  client.publish(MQTT_TOPIC_L2, String(measuredAmperL2));
  client.loop();
  client.publish(MQTT_TOPIC_L3, String(measuredAmperL3));
  
  //delay(100);
}