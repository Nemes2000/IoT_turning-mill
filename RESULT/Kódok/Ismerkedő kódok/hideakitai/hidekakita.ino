#include <Wire.h>               //built in library for iic communication
#include "MPU9250.h"            //MPU9260 sensor lib (https://github.com/hideakitai/MPU9250)

//#################################################|SETTINGS|############################################################
//----------------------------------------------|CAN BE CHANGED|---------------------------------------------------------
#define SDA 23                //SET: the pin for data transfer
#define SCK 18                //SET: the pin for clock signs
#define MPU_ADDRES 0x68       //SET: the sensor factory number

//----------------------------------------------|SHOULD NOT CHANGE|-----------------------------------------------------
MPU9250 mpu9250;              //create an mpu9250 object
MPU9250Setting setting;       //assists for the sensor set up
//######################################################################################################################

void setup() {
  Serial.begin(115200);

  //set up the iic communication
  Wire.begin(SDA,SCK);

  //settings for mpu9250
  setting.accel_fs_sel = ACCEL_FS_SEL::A2G;
  setting.fifo_sample_rate = FIFO_SAMPLE_RATE::SMPL_1000HZ;
  setting.accel_fchoice = 0x01;
  setting.accel_dlpf_cfg = ACCEL_DLPF_CFG::DLPF_45HZ;

  //set up the sensor
  mpu9250.setup(MPU_ADDRES, setting, Wire);

  //calibrate, must hold still!!
  mpu9250.calibrateAccelGyro();
}

void loop() {
  //refreshing the mpu stored data
  mpu9250.update_accel_gyro();

  Serial.println(mpu9250.getAccX());
}