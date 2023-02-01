#include <mpu9250.h>
#include <Wire.h>

bfs::Mpu9250 imu;

void setup() {
  Serial.begin(115200);
  while(!Serial) {Serial.println("Stopped");}
  
  Wire.begin(23, 18);
  //Wire.setClock(400000);
  imu.Config(&Wire, 0x68);

  if (!imu.Begin()) {
    Serial.println("Error initializing communication with IMU");
    while(1) {}
  }

  if (!imu.ConfigSrd(19)) {
    Serial.println("Error configured SRD");
    while(1) {}
  }
}


void loop() {
  if(imu.Read()){
    Serial.println(imu.accel_x_mps2());
  }
}