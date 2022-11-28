//##############################|SETTINGS|####################################################################
#define MIC_PIN 36              //SET: the pin where you get the analog input from microphon
#define ASSIST_NUM 1800         //SET: just an assist number, doesn't necessary
//############################################################################################################

void setup() {
  Serial.begin(115200);
}

void loop() {
    Serial.print("Data: ");

    //read from the pin and print that out
    Serial.print(analogRead(MIC_PIN));

    Serial.print(" , ");
    Serial.print("assist: ");
    Serial.println(ASSIST_NUM);
}
