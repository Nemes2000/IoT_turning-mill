TaskHandle_t Task1;       //creating a new thread object
TaskHandle_t Task2;       //creating a new thread object: doesn't need because by default run on core one

void setup() {
  Serial.begin(115200);

  //set up the new thread
  xTaskCreatePinnedToCore(
    loopForTemperature, /* Function to implement the task */
    "Task1",            /* Name of the task */
    10000,              /* Stack size in words */
    NULL,               /* Task input parameter */
    0,                  /* Priority of the task -> higher number higher priority*/
    &Task1,             /* Which task handle it. */
    0);                 /* Core where the task should run */
}

void loop() { 
  //main thread -> by default the thread1
  Serial.print("running on core ");

  //print out the thread id
  Serial.println(xPortGetCoreID());
  delay(5000);
}

void loopForTemperature(void * pvParameters){
  //creating an infinite loop so it will work as the main loop function
  for(;;) {
    //side thread -> by default the thread2
    Serial.print("running on core ");

    //print out the thread id
    Serial.println(xPortGetCoreID());
    delay(5000);
  }
}
