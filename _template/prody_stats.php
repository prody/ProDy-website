<?php
$servername = "localhost";
$username = "prody";
$password = "I'm Protein Dynamics";
$dbname = "prody";

// Create connection
$conn = new mysqli($servername, $username, $password, $dbname);
// Check connection
if ($conn->connect_error) {
    die("Connection failed: " . $conn->connect_error);
}

$sql = "SELECT SUM(number) AS value_sum FROM downloads";
$result = $conn->query($sql);

$row = $result->fetch_assoc();
$sum = $row["value_sum"];

$result->free();
$conn->close();

echo number_format($sum);
?>
