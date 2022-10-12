const int micPin = 36;
const int mid = 1800;

void setup() {
  Serial.begin(9600);
}

void loop() {
    Serial.print("Jel: ");
    Serial.print(analogRead(micPin));
    Serial.print(" , ");
    Serial.print("allando: ");
    Serial.println(mid);
}
