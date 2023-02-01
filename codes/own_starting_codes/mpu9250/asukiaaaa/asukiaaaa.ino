#include <Wire.h>                   //built in library for iic communication
#include <MPU9250_asukiaaa.h>       //MPU9260 sensor lib (https://github.com/hideakitai/MPU9250)

//############################################|SETTINGS|########################################################################################
//-----------------------------------------|Can be setted|--------------------------------------------------------------------------------------
#define SDA_PIN 23                    //SET: the pin number where the iic sda is located
#define SCK_PIN 18                    //SET: the pin number where the iic sck is located

//---------------------------------------|Should not change|------------------------------------------------------------------------------------
MPU9250_asukiaaa mpu9250;              //create an mpu9250 object
//##############################################################################################################################################

void setup() {
  Serial.begin(115200);

  //set up the iic communication
  Wire.begin(SDA, SCK);

  //set up the sensor
  mpu9250.setWire(&Wire);

  //configure accelometer
  mpu9250.beginAccel();
}

void loop() {
  //refresh the sensor's accelometer datas
  mpu9250.accelUpdate();

  Serial.print(mpu9250.accelX());
  Serial.print(", ");
  Serial.print(mpu9250.accelY());
  Serial.print(", ");
  Serial.print(mpu9250.accelZ());
  Serial.print(", ");
  Serial.print(mpu9250.accelSqrt());
  Serial.println();
}









