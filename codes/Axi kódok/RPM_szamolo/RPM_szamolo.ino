#include "EspMQTTClient.h"    // MQTT lib (https://github.com/plapointe6/EspMQTTClient)

// ###################################|MQTT KAPCSOLAT KIALAKITASA|#####################################################
// WiFi beallitasok
EspMQTTClient client(
  "",                   // BSID
  "",                   // AC Jelszo
  "",                   // MQTT Szerver IP cime
  "",                   // MQTT Szerver Felhasznaloven (Ha nem kell a szerernek felhasznalao nevet megadni akkor TOROLD KI EZT A SORT)
  "",                   // MQTT Szerver Jelszava (Ha nem kell a szerernek jelszot megadni akkor TOROLD KI EZT A SORT)
  "ESP-12E-RPM",        // MQTT kliens neve (ESP neve, tetszoleges)
  1883                  // MQTT szerver portja
);
// ####################################################################################################################


// ########################################|MERESI BEALLITASOK|########################################################
#define OPTO_PIN D5           // Opto pinje
#define TIMEOUT 0.4           // Hany masodperc semmitteves utan timeout-oljon
#define FINOM_EJTES true      // Ha timeout-ol a fordulatszam mero akkor fokozatosan eljtse-e a jelet vagy sem
#define FINOM_SZORZO 0.4      // Ha finom ejtest alkalmazunk, akkor milyen gyorsan csengjen le a jel
#define MQTT_TOPIC "ESP/RPM"  // Az a csatorna/topic ahova kuldi az adatot

long elotte, most;            // Az elozo es mostani interrupt ideje microsec-ben
int tombIndex = 0;            // a lyukakKozottiIdo indexe, ora modszerrel (opre) toltjuk fel/frissitjuk az ertekeket a tombben
long tombOsszeg = 0.0;        // Ebbe szamoljuk mindig a tombnek az osszeget
double RPM = 0.0;             // Fordulatszam
bool timeoutVan = false;      // Van-e timeout?
long lyukakKozottiIdo[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};    // Ablakkos modszer miatt, eltaroljuk a lyukak kozotti erteket es az alapjan szamolunk
                                                              // Igy a korulfordulasi ido ennek a tomb elemeinek az osszege
// ####################################################################################################################

// ##########################################|SEGED FUGGVENYEK|########################################################
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

void onConnectionEstablished()
{
  //client.subscribe(MQTT_TOPIC, [](const String & payload) {
  //  Serial.println(payload);
  //});
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

  // Init vege
  Serial.println("Initialization Finished!");
  
  elotte = micros();
}

// ########################################|LOOP|######################################################################
void loop() {

  // MQTT-hez
  client.loop();

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
  client.publish(MQTT_TOPIC, String(RPM));
}
