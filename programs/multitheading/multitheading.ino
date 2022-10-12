

TaskHandle_t Task1;
//TaskHandle_t Task2; //dont need because by default run on core one

void setup() {
  Serial.begin(115200);

  xTaskCreatePinnedToCore(
    loopForTemperature, /* Function to implement the task */
    "Task1",            /* Name of the task */
    10000,              /* Stack size in words */
    NULL,               /* Task input parameter */
    0,                  /* Priority of the task -> higher number higher priority*/
    &Task1,             /* Task handle. */
    0);                 /* Core where the task should run */

  //xTaskCreatePinnedToCore( loopForAccelerometer, "Task2", 10000, NULL, 0, &Task2, 1); 
}

void loop() { //for accelerometer
  Serial.print("running on core ");
  Serial.println(xPortGetCoreID());
  delay(5000);
}

void loopForTemperature(void * pvParameters){
  for(;;) { //infinete loop -> like loop on core 1
    Serial.print("running on core ");
    Serial.println(xPortGetCoreID());
    delay(5000);
  }
}
