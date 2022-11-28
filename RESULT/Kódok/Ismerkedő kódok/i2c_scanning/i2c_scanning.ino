#include <Wire.h>        //built in iic library

//##########################|SETTINGS|#######################################################
#define SDA 23        //SET: the pin where you get data
#define SCK 18        //SET: the pin where you send the clock signs
//###########################################################################################

void setup()
{
  Serial.begin(9600);

  //set up the iic communication
  Wire.begin(SDA,SCK);

  //wait for open the serila monnitor
  while (!Serial);   

  Serial.println("\nI2C Scanner");
}
 
void loop()
{
  byte error, address;
  int nDevices;
 
  Serial.println("Scanning...");
 
  nDevices = 0;

  //go through all the possible addresses (the address is a byte long)
  for(address = 1; address < 127; address++ )
  {
    // The i2c_scanner uses the return value of
    // the Write.endTransmisstion to see if
    // a device did acknowledge to the address.
    Wire.beginTransmission(address);
    error = Wire.endTransmission();
 
    if (error == 0)
    {
      Serial.print("I2C device found at address 0x");
      if (address<16)
        Serial.print("0");
      Serial.print(address,HEX);
      Serial.println("  !");
 
      nDevices++;
    }
    else if( error == 1){
      Serial.println("data exceeds buffer capacity");
    }
    else if( error == 2){
      Serial.print("address byte followed by NACK");
    }
    else if(error == 3){
      Serial.println("data byte followed by NACK");
    }
    else if (error==4)
    {
      Serial.print("Unknown error at address 0x");
      if (address<16)
        Serial.print("0");
      Serial.println(address,HEX);
    }  
    else if(error == 5){
      Serial.println("timeout");
    }  
  }
  if (nDevices == 0)
    Serial.println("No I2C devices found\n");
  else
    Serial.println("done\n");
 
  delay(5000);           // wait 5 seconds for next scan
}